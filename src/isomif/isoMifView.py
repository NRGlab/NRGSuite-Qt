import sys
import os

probesLab = ["HYD", "ARM", "DON", "ACC", "POS", "NEG"]
matchIn = ""
outDir = "./"
prefix = ""
cg = 0
res = 1
tcg = 0
ca = []
pseudo = []
rot = []
cen = []
va = []
vb = []
mifV1int = []
mifV2int = []
sm = "taninorm"
p1Path = None
p2Path = None
m1Path = None
m2Path = None

# Parse command line arguments
i = 0
while i < len(sys.argv):
    if sys.argv[i] == "-m":
        matchIn = sys.argv[i+1]
    elif sys.argv[i] == "-o":
        outDir = sys.argv[i+1]
    elif sys.argv[i] == "-p":
        prefix = sys.argv[i+1]
    elif sys.argv[i] == "-p1":
        p1Path = sys.argv[i+1]
    elif sys.argv[i] == "-p2":
        p2Path = sys.argv[i+1]
    elif sys.argv[i] == "-m1":
        m1Path = sys.argv[i+1]
    elif sys.argv[i] == "-m2":
        m2Path = sys.argv[i+1]
    elif sys.argv[i] == "-g":
        cg = int(sys.argv[i+1])
    elif sys.argv[i] == "-s":
        sm = sys.argv[i+1]
    elif sys.argv[i] == "-h":
        print("##################\nWelcome to pipeIsoMifView\n##################")
        print("-m         <path to isoMif file>")
        print("-o         <isoMifView output directory>")
        print("-p         <prefix for output files>")
        print("-p1        <protein 1 path>")
        print("-p2        <protein 2 path>")
        print("-m1        <mif 1 path>")
        print("-m2        <mif 2 path>")
        print("-g         <coarse grain step>")
        print("-s         <similarity measure to find best clique>")
        print("-h         <print help menu>")
        sys.exit()
    i += 1

res = cg if cg > 0 else 1

if outDir == "":
    exit('no output_dir specified with -o')

probeNames = []
pbColors = ["aquamarine", "brightorange", "blue", "red", "limegreen", "lightmagenta"]

best = {}
cc = 0

# Parse the match file if similarity measure is given
if sm:
    with open(matchIn, 'r') as f:
        for line in f:
            if line.startswith("REMARK CLIQUE CG"):
                parts = line.split()
                if parts[3] == str(cg):
                    data = {
                        "nodes": int(parts[5]),
                        "nodesm": int(parts[7]),
                        "nodesmw": float(parts[9]),
                        "normnodes": float(parts[11]),
                        "normnodesrmsd": float(parts[13]),
                        "tani": float(parts[15]),
                        "tanim": float(parts[17]),
                        "tanimw": float(parts[19]),
                        "taninorm": float(parts[21]),
                        "nrg": float(parts[23]),
                        "ss1": int(parts[25]),
                        "ss2": int(parts[27]),
                        "ss1m": int(parts[29]),
                        "ss2m": int(parts[31]),
                        "ligrmsd": float(parts[33])
                    }

                    if cc == 0:
                        for key in data:
                            best[key] = [cc, data[key]]
                    else:
                        for key, value in data.items():
                            if (key in ["nodes", "nodesm", "normnodes", "normnodesrmsd", "tani", "tanim", "taninorm", "nrg"]
                                    and value > best[key][1]) or (key == "ligrmsd" and value < best[key][1]):
                                best[key] = [cc, value]

                    cc += 1

    print(f"Best {sm}: {best[sm][0]} {best[sm][1]}")

# Read the match file to retrieve nodes and other information
cc = 0
flagStore = False

def get_data(matchIn):
    with open(matchIn, 'r') as f:
        max_pb = 0
        for line in f:
            line = line.strip()

            if "REMARK CLIQUE CG" in line:
                temp_tcg = int(line.split()[3])

            if not line.startswith("REMARK"):
                pb = int(line.split()[0])
                if pb > max_pb:
                    max_pb = pb
    max_pb = 5
    return temp_tcg, max_pb

temp_tcg, max_pb = get_data(matchIn)

data = [["" for _ in range(max_pb+1)] for _ in range(temp_tcg+1)]

with open(matchIn, 'r') as f:
    for line in f:
        line = line.strip()
        if "mif_file_1:" in line:
            if not m1Path:
                m1Path = line.split(":")[1].strip()
            if not p1Path:
                p1Path = m1Path.replace(".mif", "_cpy.pdb")
            mif1 = os.path.basename(m1Path)

        if "mif_file_2:" in line:
            if not m2Path:
                m2Path = line.split(":")[1].strip()
            if not p2Path:
                p2Path = m2Path.replace(".mif", "_cpy.pdb")
            mif2 = os.path.basename(m2Path)

        if "REMARK CLIQUE CG" in line:
            tcg = int(line.split()[3])
            if flagStore:
                break
            flagStore = sm and best[sm][0] == cc
            cc += 1 if tcg == cg else 0

        if "REMARK ROTMAT" in line and tcg == cg:
            rot = [list(map(float, line.split()[2:5])),
                   list(map(float, line.split()[5:8])),
                   list(map(float, line.split()[8:11]))]

        if "REMARK CENTRES" in line and tcg == cg:
            cen = [list(map(float, line.split()[2:5])),
                   list(map(float, line.split()[5:8]))]

        if not line.startswith("REMARK") and tcg == cg:
            if flagStore == 1:
                l = line.split()  # Split the line into a list based on whitespace
                pb = int(l[0])  # Get the first element as pb
                if data[tcg][pb] == "":
                    data[tcg][pb] = f"{l[1]};{l[2]};{l[3]};{l[4]};{l[5]};{l[6]}\n"
                else:
                    data[tcg][pb] += f"{l[1]};{l[2]};{l[3]};{l[4]};{l[5]};{l[6]}\n"

# Additional processing functions for storing and manipulating MIF data
def store_mif(mif_file, mif_data):
    with open(mif_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#probe"):
                probe_name = line.split()[-1]
                probeNames.append(probe_name)
            elif line.startswith("#") or not line:
                continue
            else:
                info = list(map(float, line.split()))
                for i in range(3, len(info) - 5, 3):
                    pbid = (i // 3) - 1
                    if int(info[i]) == 1:
                        temp_list = info[:3] + info[-5:-1]
                        for item in temp_list:
                            mif_data[pbid].append(item)

mifV1 = [[] for _ in range(len(probesLab))]
mifV2 = [[] for _ in range(len(probesLab))]
if m1Path and os.path.exists(m1Path):
    store_mif(m1Path, mifV1)
if m2Path and os.path.exists(m2Path):
    store_mif(m2Path, mifV2)

def print_mif(mif_id, mif_data, mif_int):
    it = 0
    pdb_str = "cmd.read_pdbstr(\"\"\""
    mif_int = [[[[] for _ in range(2)] for _ in range(res+1)] for _ in range(len(mif_data))]
    for i, probe in enumerate(mif_data):
        if len(probe):
            mif_int[i][res][0] = it
            for j in range(0, len(probe), 7):
                if probe[j + 3 + res] == 1:
                    coor = probe[j:j+3]
                    ncoor = coor.copy()
                    if mif_id == 1:
                        for k in range(3):
                            ncoor[k] = cen[1][k]
                            for l in range(3):
                                ncoor[k] += (coor[l] - cen[0][l]) * rot[k][l]
                    pdb_str += f"HETATM{it:5d}  N   {i:3} A0000    {ncoor[0]:8.3f}{ncoor[1]:8.3f}{ncoor[2]:8.3f}  0.00 10.00           N\\\n"
                    it += 1
            mif_int[i][res][1] = it - 1
    pdb_str += "TER \\\n\"\"\",\"{}_mif{}\")\n".format(tag, mif_id)
    return pdb_str, mif_int

mif1 = mif1.replace(".mif", "")
mif2 = mif2.replace(".mif", "")
tag = f"{mif1}_{mif2}"
if prefix:
    tag = f"{prefix}_{tag}"

mif1str, mifV1int = print_mif(1, mifV1, mifV1int)
mif2str, mifV2int = print_mif(2, mifV2, mifV2int)

# Create protein 1 file
p1str = "cmd.read_pdbstr(\"\"\""
with open(p1Path, 'r') as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            coor = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
            ncoor = [cen[1][i] + sum((coor[j] - cen[0][j]) * rot[i][j] for j in range(3)) for i in range(3)]
            p1str += f"{line[:30]}{ncoor[0]:8.3f}{ncoor[1]:8.3f}{ncoor[2]:8.3f}{line[54:-1]}\\\n"
p1str += f"TER \\\n\"\"\",\"{mif1}\")\n"

p2str = 'cmd.read_pdbstr("""'
with open(p2Path, 'r') as f:
    for line_counter, line in enumerate(f):
        line = line.strip()  # Equivalent to chomp in Perl
        if line_counter != 0 and line != 'TER' and not line.startswith('CONECT'):
            while len(line) < 80:
                line += ' '
        if line == 'TER':
            line += '   '
        p2str += line + "\\\n"
p2str += 'TER \\\n""","' + mif2 + '")\n'

def print_mif_pml(mif_v, mif_v_int, probe_num, tag, sphere_scale, full_tag):
    pml_str = ""
    # Loop through the items in mif_v
    for i in range(len(mif_v)):
        if mif_v[i]:  # Equivalent to checking if @mifV1[i] is true
            # Assuming res is defined or passed as a parameter.
            if mif_v_int[i][res][0] != mif_v_int[i][res][1]:
                pml_str += (f"create mif_{tag}_{probeNames[i]}, id {mif_v_int[i][res][0]}-{mif_v_int[i][res][1]} & "
                            f"{full_tag}_mif{probe_num}\n")
                pml_str += (f"show spheres, mif_{tag}_{probeNames[i]}\n"
                            f"set sphere_scale,{sphere_scale},mif_{tag}_{probeNames[i]}\n"
                            f"set sphere_transparency,0.6,mif_{tag}_{probeNames[i]}\n"
                            "rebuild\n")
                pml_str += f"color {pbColors[i]},mif_{tag}_{probeNames[i]}\n"
                pml_str += f"hide nonbonded,mif_{tag}_{probeNames[i]}\n"

    pml_str += f"delete {full_tag}_mif{probe_num}\n"
    return pml_str

# Define NPML output file
with open(f"{outDir}{tag}.pml", "w") as NPML:
    # Initial PML settings
    NPML.write(f"{p1str}{p2str}show cartoon\n")
    NPML.write(f"hide lines\n")
    NPML.write(f"set connect_mode,1\n")
    NPML.write(f"{mif1str}{mif2str}")

    if cg == -1:
        for nod in ca:
            s = nod.split(";")
            NPML.write(f"show lines, resi {s[1]} & chain {s[2]} & {mif1}\n")
            NPML.write(f"show lines, resi {s[4]} & chain {s[5]} & {mif2}\n")
    elif cg == -3:
        ps1 = ""
        ps2 = ""
        for p in pseudo:
            s = p.split(";")
            id = 0
            coor = [s[1], s[2], s[3]]
            ncoor = [0, 0, 0]
            for i in range(3):
                ncoor[i] = cen[1][i]
                for k in range(3):
                    ncoor[i] += (float(coor[k]) - cen[0][k]) * rot[i][k]

            ps1 += f"HETATM{id:5d}  CA  {s[0].upper()} A        {ncoor[0]:8.3f}{ncoor[1]:8.3f}{ncoor[2]:8.3f}  0.00 10.00           C  \\\n"
            ps2 += f"HETATM{id:5d}  CA  {s[4].upper()} A        {s[5]:8.3f}{s[6]:8.3f}{s[7]:8.3f}  0.00 10.00           C  \\\n"
            id += 1

        NPML.write(f'cmd.read_pdbstr("""{ps1}TER \\\n""", "{tag}_1_pseudo")\n')
        NPML.write(f'cmd.read_pdbstr("""{ps2}TER \\\n""", "{tag}_2_pseudo")\n')
    elif cg == -2:
        id = 0
        NPML.write('set connect_mode,1\ncmd.read_pdbstr("""')

        # First set of nodes
        for p in range(6):
            if len(va[p]) > 0:
                start = id
                for i in range(0, len(va[p]), 3):
                    NPML.write(
                        f"HETATM{id:5d}  CA  NRG A        {va[p][i]:8.3f}{va[p][i + 1]:8.3f}{va[p][i + 2]:8.3f}  0.00 10.00           C  \\\n")
                    id += 1
                stop = id - 1
                NPML.write(f'create {probeNames[p]}_{mif1}, id {start}-{stop} & {tag}_1_nodes\n')
                NPML.write(f'set sphere_scale,0.25,{probeNames[p]}_{mif1}\n')
                NPML.write(f'show spheres, {probeNames[p]}_{mif1}\nrebuild\n')
                NPML.write(f'color {pbColors[p]}, {probeNames[p]}_{mif1}\n')

        NPML.write(f'TER \\\n""", "{tag}_1_nodes")\n')

        # Second set of nodes
        id = 0
        NPML.write('set connect_mode,1\ncmd.read_pdbstr("""')
        for p in range(6):
            if len(vb[p]) > 0:
                start = id
                for i in range(0, len(vb[p]), 3):
                    NPML.write(
                        f"HETATM{id:5d}  CA  NRG A        {vb[p][i]:8.3f}{vb[p][i + 1]:8.3f}{vb[p][i + 2]:8.3f}  0.00 10.00           C  \\\n")
                    id += 1
                stop = id - 1
                NPML.write(f'create {probeNames[p]}_{mif2}, id {start}-{stop} & {tag}_2_nodes\n')
                NPML.write(f'set sphere_scale,0.15,{probeNames[p]}_{mif2}\n')
                NPML.write(f'show spheres, {probeNames[p]}_{mif2}\nrebuild\n')
                NPML.write(f'color {pbColors[p]}, {probeNames[p]}_{mif2}\n')

        NPML.write('TER \\\n""", "{tag}_2_nodes")\n')
    else:
        #TODO: previous cases not tested
        ids = 0
        str1 = 'cmd.read_pdbstr("""'
        str2 = 'cmd.read_pdbstr("""'
        strSel = ""
        for j in range(6):  # For each probe
            nodes = data[cg][j].split("\n")[:-1]
            if nodes:
                start = ids
                for node in nodes:
                    info = list(map(float, node.split(";")))
                    coor = [float(info[0]), float(info[1]), float(info[2])]
                    ncoor = [cen[1][i] + sum((coor[k] - cen[0][k]) * rot[i][k] for k in range(3)) for i in range(3)]

                    str1 += f"HETATM{ids:5d}  CA  NRG A        {ncoor[0]:8.3f}{ncoor[1]:8.3f}{ncoor[2]:8.3f}  0.00 10.00           C  \\\n"
                    str2 += f"HETATM{ids:5d}  CA  NRG A        {info[3]:8.3f}{info[4]:8.3f}{info[5]:8.3f}  0.00 10.00           C  \\\n"
                    ids += 1

                stop = ids - 1
                strSel += f'create {probeNames[j]}_{mif1}, id {start}-{stop} & {tag}_1_nodes\n'
                strSel += f'set sphere_scale,0.25,{probeNames[j]}_{mif1}\nshow spheres, {probeNames[j]}_{mif1}\nrebuild\n'
                strSel += f'color {pbColors[j]}, {probeNames[j]}_{mif1}\n'
                strSel += f'create {probeNames[j]}_{mif2}, id {start}-{stop} & {tag}_2_nodes\n'
                strSel += f'set sphere_scale,0.15,{probeNames[j]}_{mif2}\nshow spheres, {probeNames[j]}_{mif2}\nrebuild\n'
                strSel += f'color {pbColors[j]}, {probeNames[j]}_{mif2}\n'

        str1 += f'TER \\\n""","{tag}_1_nodes")\n'
        str2 += f'TER \\\n""","{tag}_2_nodes")\n'
        NPML.write(str1 + str2 + strSel)

    mstr1 = print_mif_pml(mifV1, mifV1int, 1, mif1, 0.25, tag)
    mstr2 = print_mif_pml(mifV2, mifV2int, 2, mif2, 0.15, tag)

    NPML.write(f"set connect_mode,1\n{mstr1}{mstr2}")

    NPML.write(f"remove hydrogens\nshow sticks, HET\ndelete {tag}_1_nodes\ndelete {tag}_2_nodes\n")

    pseudo_lab = ["HYD", "ARM", "DON", "ACC", "DOA"]

    NPML.write(f"set sphere_scale, 0.25, {tag}_1_pseudo\n"
               f"color aquamarine, resn HYD & {tag}_1_pseudo\n"
               f"color brightorange, resn ARM & {tag}_1_pseudo\n"
               f"color blue, resn DON & {tag}_1_pseudo\n"
               f"color red, resn ACC & {tag}_1_pseudo\n"
               f"color limegreen, resn DOA & {tag}_1_pseudo\n"
               f"show spheres, {tag}_1_pseudo\n"
               f"set sphere_scale, 0.15, {tag}_2_pseudo\n"
               f"color aquamarine, resn HYD & {tag}_2_pseudo\n"
               f"color brightorange, resn ARM & {tag}_2_pseudo\n"
               f"color blue, resn DON & {tag}_2_pseudo\n"
               f"color red, resn ACC & {tag}_2_pseudo\n"
               f"color limegreen, resn DOA & {tag}_2_pseudo\n"
               f"show spheres, {tag}_2_pseudo\n")

    # Assuming you have opened NPML file, you will close it at the end of the script
NPML.close()
