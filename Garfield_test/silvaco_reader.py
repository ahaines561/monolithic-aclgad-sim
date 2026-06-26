import sys

QUANTITY_NAMES = {
    513: ("Region", ""),
    5:   ("Boron", "cm^-3"),
    151: ("ActiveBoron", "cm^-3"),
    3:   ("Phosphorus", "cm^-3"),
    152: ("ActivePhosphorus", "cm^-3"),
    71:  ("DonorConcentration", "cm^-3"),
    72:  ("AcceptorConcentration", "cm^-3"),
    115: ("NetDoping", "cm^-3"),
    149: ("TotalDoping", "cm^-3"),
    100: ("ElectrostaticPotential", "V"),
    103: ("ElectricFieldMagnitude", "V/cm"),
    120: ("ElectricField_x", "V/cm"),
    121: ("ElectricField_y", "V/cm"),
    124: ("LatticeTemperature", "K"),
}
UM_TO_CM = 1.0e-4

class SilvacoMesh:
    def __init__(self):
        self.version = None
        self.counts = { } # nCoord, nNode, nEdge, nReg, nTri
        self.bbox = {} # xmin, xmax, ymin, ymax
        self.k_is_mesh_header = None  # True for ATHENA .str, False for ATLAS .sta
        self.k_raw = None # raw 'k' fields when not a mesh header
        self.vertices = {}  # id -> (x, y, z) in cm
        self.edges = []  # (id, v1, v2, flag)
        self.triangles = []  # (id, region, (v1,v2,v3), (n1,n2,n3))
        self.regions = {} # region_id -> material_code
        self.material_table = [] # from 'j' record
        self.codes = None  # solution quantity codes (column order)
        self.node_vals = {} # node_id -> [values] matching self.codes

def parse(path, mesh=None):
    if mesh is None:
        mesh = SilvacoMesh()
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            tag = line[0]
            if tag == "v":
                mesh.version = line[1:].strip()
            elif tag == "j":
                p = line.split()
                n = int(p[1])
                mesh.material_table = [int(x) for x in p[2:2 + n]]
            elif tag == "k":
                p = line.split()
                try:
                    mesh.counts = dict(
                        nCoord=int(p[1]), nNode=int(p[2]), nEdge=int(p[3]),
                        nReg=int(p[4]), nTri=int(p[5]),
                    )
                    mesh.bbox = dict(
                        xmin=0.0, xmax=float(p[7]),
                        ymin=float(p[8]), ymax=float(p[9]),
                    )
                    mesh.k_is_mesh_header = True
                except (IndexError, ValueError):
                    mesh.k_raw = p[1:]
                    mesh.k_is_mesh_header = False

            elif tag == "c":
                p = line.split()
                cid = int(p[1])
                x = float(p[2]) * UM_TO_CM
                y = float(p[3]) * UM_TO_CM
                z = float(p[4]) * UM_TO_CM if len(p) > 4 else 0.0
                mesh.vertices[cid] = (x, y, z)

            elif tag == "e" and line[1] == " ":
                p = line.split()
                mesh.edges.append(
                    (int(p[1]), int(p[2]), int(p[3]),
                     int(p[4]) if len(p) > 4 else 0))

            elif tag == "t" and line[1] == " ":
                p = line.split()
                tid = int(p[1])
                region = int(p[2])
                verts = (int(p[3]), int(p[4]), int(p[5]))
                neigh = (int(p[6]), int(p[7]), int(p[8])) if len(p) >= 9 \
                    else (0, 0, 0)
                mesh.triangles.append((tid, region, verts, neigh))

            elif tag == "r" and line[1] == " ":
                p = line.split()
                mesh.regions[int(p[1])] = int(p[2])

            elif tag == "s" and line[1] == " ":
                p = line.split()
                ncode = int(p[1])
                mesh.codes = [int(x) for x in p[2:2 + ncode]]

            elif tag == "n" and mesh.codes is not None and line[1] == " ":
                p = line.split()
                field_a = int(p[2])
                if field_a == 0:
                    continue
                nid = int(p[1])
                raw = p[3:3 + len(mesh.codes)]
                if len(raw) < len(mesh.codes):
                    continue
                try:
                    vals = [float(x) for x in raw]
                except ValueError:
                    continue
                mesh.node_vals[nid] = vals
    return mesh

def summarize(mesh):
    print("Silvaco mesh summary")
    print("  generator     :", mesh.version,
          "(.str / ATHENA mesh header)" if mesh.k_is_mesh_header
          else "(.sta / ATLAS solution header)")
    if mesh.k_is_mesh_header:
        print("  declared counts:", mesh.counts)
        print("  bounding box  : x=[{xmin}, {xmax}]  y=[{ymin}, {ymax}] (um)"
              .format(**{k: round(v, 3) for k, v in mesh.bbox.items()}))
    else:
        print("  k header (raw):", mesh.k_raw,
              "  <- solution metadata, not mesh counts")
    print("  vertices read :", len(mesh.vertices))
    print("  triangles read:", len(mesh.triangles))
    print("  edges read    :", len(mesh.edges))
    print("  regions       :", len(mesh.regions),
          "->", dict(sorted(mesh.regions.items())))
    if mesh.codes:
        print("  solution columns ({} quantities):".format(len(mesh.codes)))
        for i, c in enumerate(mesh.codes):
            name, units = QUANTITY_NAMES.get(c, ("code%d" % c, ""))
            mark = "  <-- FIELD" if c in (120, 121) else ""
            print("      col {:2d} : code {:4d} : {} [{}]{}"
                  .format(i, c, name, units, mark))
        print("  nodal records :", len(mesh.node_vals))
        # for athena
    if mesh.k_is_mesh_header:
        print("\nconsistency vs. 'k' header:")
        chk = [
            ("vertices", len(mesh.vertices), mesh.counts.get("nCoord")),
            ("triangles", len(mesh.triangles), mesh.counts.get("nTri")),
            ("edges", len(mesh.edges), mesh.counts.get("nEdge")),
            ("regions", len(mesh.regions), mesh.counts.get("nReg")),
        ]
        for name, got, declared in chk:
            flag = "OK" if got == declared else "MISMATCH"
            print("  {:10s} read {:>8} declared {:>8}  [{}]"
                  .format(name, got, declared, flag))

def sample_triangle(mesh):
    if not mesh.triangles:
        return
    tid, region, verts, neigh = mesh.triangles[0]
    print("\nexample element t{} (region {}):".format(tid, region))
    for v in verts:
        if v in mesh.vertices:
            x, y, z = mesh.vertices[v]
            print("    vertex {:>6} at ({:.4f}, {:.4f}) um"
                  .format(v, x / UM_TO_CM, y / UM_TO_CM))
    print("    neighbour triangles:", neigh,
          "(negative = boundary)")
def spot_check_values(mesh):
    if not mesh.codes or not mesh.node_vals:
        return
    col = {c: i for i, c in enumerate(mesh.codes)}

    def show(node_id):
        if node_id not in mesh.node_vals:
            return
        vals = mesh.node_vals[node_id]
        cid = node_id + 1 
        xyz = mesh.vertices.get(cid)
        loc = ("({:.2f}, {:.2f}) um".format(xyz[0] / UM_TO_CM, xyz[1] / UM_TO_CM)
               if xyz else "?")
        print("  node {:>6} at {}".format(node_id, loc))
        for code, name, fmt in [
            (513, "Region", "{:.0f}"),
            (115, "NetDoping", "{:.3e}"),
            (100, "Potential", "{:.4g} V"),
            (103, "|E|", "{:.4g} V/cm"),
            (120, "Ex", "{:.4g} V/cm"),
            (121, "Ey", "{:.4g} V/cm"),
        ]:
            if code in col:
                print("      {:10s} = ".format(name)
                      + fmt.format(vals[col[code]]))

    print("\nvalue spot-check (verify columns hold sane physical numbers):")
    for nid in (0, 30000, 59684):
        show(nid)
    if 121 in col:
        ey_col = col[121]
        eys = [v[ey_col] for v in mesh.node_vals.values()]
        print("  Ey over all nodes: min={:.4g}  max={:.4g} V/cm"
              .format(min(eys), max(eys)))

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    mesh = parse(sys.argv[1])
    if len(sys.argv) > 2:
        parse(sys.argv[2], mesh)
    summarize(mesh)
    sample_triangle(mesh)
    spot_check_values(mesh)

if __name__ == "__main__":
    main()