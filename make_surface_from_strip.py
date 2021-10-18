#!/usr/bin/env python3
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import h5py
from utils import Params
from make_xdmf import header, grid_begin, grid_end, mesh_edge,\
    attrib_scalar_node, footer, timestamp, attrib_vector_node, mesh_face_quad

parser = argparse.ArgumentParser(description="Make xdmf surface from strip")
parser.add_argument("folder", type=str, help="Folder")
parser.add_argument("-t_min", type=float, default=0.0, help="t_min")
parser.add_argument("-t_max", type=float, default=np.inf, help="t_max")
parser.add_argument("-skip", type=int, default=0, help="Skip timesteps")
args = parser.parse_args()

params = Params(args.folder)
t0 = params.get_tmin()

files = os.listdir(args.folder)

posf = dict()
for file in files:
    if file[:11] == "data_from_t" and file[-3:] == ".h5":
        t = float(file[11:-3])
        posft = os.path.join(args.folder, file)
        try:
            with h5py.File(posft, "r") as h5f:
                for grp in h5f:
                    posf[float(grp)] = (posft, grp)
        except:
            pass

ts = []
for t in list(sorted(posf.keys())[::(args.skip+1)]):
    if t >= args.t_min and t <= args.t_max:
        ts.append(t)

print("Found timesteps")

possible_fields = [["u", "Vector", "Node"],
                   ["c", "Scalar", "Node"],
                   ["tau", "Scalar", "Node"],
                   ["p", "Scalar", "Node"],
                   ["rho", "Scalar", "Node"],
                   ["H", "Scalar", "Node"],
                   ["n", "Vector", "Node"] #,
                   #["dA", "Scalar", "Face"],
                   #["dA0", "Scalar", "Face"],
                   #["dl", "Scalar", "Edge"],
                   #["dl0", "Scalar", "Edge"]
                   ]

nodes_ = []
edges_ = []
# faces_ = []
data_ = dict()
time_ = []

fields = []
posft0, grp0 = posf[ts[0]]
with h5py.File(posft0, "r") as h5f:
    for field in possible_fields:
        if field[0] in h5f[grp0]:
            fields.append(field)
for field, vtype, vloc in fields:
    data_[field] = []

node_count = 0
node_count_ = [node_count]
for it, t in enumerate(ts):
    posft, grp = posf[t]
    with h5py.File(posft, "r") as h5f:
        nodes = np.array(h5f[grp + "/points"])
        has_faces = grp + "/faces" in h5f
        if has_faces:
            exit("No support for faces")
        edges = []
        has_edges = grp + "/edges" in h5f
        if has_edges:
            edges = np.array(h5f[grp + "/edges"])

        nodes_.append(nodes)
        edges_.append(edges + node_count)
        node_count += len(nodes)
        node_count_.append(node_count)
        time_.append(t * np.ones((len(nodes), 1)))
        for field, vtype, vloc in fields:
            data_loc = np.array(h5f[grp + "/" + field])
            data_[field].append(data_loc)

node2next = dict()
for c1_, c2_, nc1, nc2 in zip(data_["c"][:-1], data_["c"][1:], node_count_[:-1], node_count_[1:]):
    for i1, c1 in enumerate(c1_):
        i2 = next((k for k, l in enumerate(c2_[i1:]) if l==c1), None)
        # print(i1+nc1, i2+nc2)
        if i2 is not None:
            node2next[i1+nc1] = i1+i2+nc2

#along_edges = np.zeros((node_count_[-2], 2), dtype=int)
#for i, (v1, v2) in enumerate(node2next.items()):
#    along_edges[i, :] = [v1, v2]
#print(along_edges)

nodes = np.vstack(nodes_)
# edges = np.vstack(edges_ + [along_edges])
edges = np.vstack(edges_)
time = np.vstack(time_)
data = dict()
for field, vtype, vloc in fields:
    data[field] = np.vstack(data_[field])
c = data["c"]

node2r = dict()
node2l = dict()
for edge in edges:
    i1, i2 = edge
    imax, imin = (i1, i2) if c[i1] > c[i2] else (i2, i1)
    node2l[imax] = imin
    node2r[imin] = imax

faces = []
for i0 in range(len(nodes)):
    if i0 in node2next and i0 in node2r and node2r[i0] in node2next:
        i1 = node2r[i0]
        i2 = node2next[i1]
        i3 = node2next[i0]
        faces.append([i0, i1, i2, i3])
faces = np.array(faces, dtype=int)
#print(faces)
#exit()

#print(nodes)
#print(edges)
#print(time)
for field, vtype, vloc in fields:
    print(field, data[field].shape)
print("time", time.shape)

h5filename = "time_data.h5"
with h5py.File(os.path.join(args.folder, h5filename), "w") as h5f:
    h5f.create_dataset("surface/points", data=nodes, dtype='float64')
    # h5f.create_dataset("surface/edges", data=edges, dtype='int64')
    h5f.create_dataset("surface/faces", data=faces, dtype='int64')
    h5f.create_dataset("surface/time", data=time, dtype='float64')
    for field, vtype, vloc in fields:
        h5f.create_dataset("surface/" + field, data=data[field], dtype='float64')

text = header
text += grid_begin
#text += mesh_edge.format(num_edges=len(edges), num_nodes=len(nodes),
#                         filename=h5filename, edges_loc="surface/edges",
#                         nodes_loc="surface/points")
text += mesh_face_quad.format(num_faces=len(faces), num_nodes=len(nodes),
                         filename=h5filename, faces_loc="surface/faces",
                         nodes_loc="surface/points")
text += timestamp.format(time=0.)
attrib = ""
attrib += attrib_scalar_node
text += attrib.format(num_nodes=len(nodes),
                      # num_faces=len(faces),
                      num_edges=len(edges),
                      filename=h5filename,
                      field_loc="surface/time",  # +field,
                      field="time")  #field)
for field, vtype, vloc in fields:
    attrib = ""
    if vtype == "Scalar":
        attrib += attrib_scalar_node
    elif vtype == "Vector":
        attrib += attrib_vector_node
    text += attrib.format(num_nodes=len(nodes),
                          # num_faces=len(faces),
                          num_edges=len(edges),
                          filename=h5filename,
                          field_loc="surface/" + field,
                          field=field)

text += grid_end
text += footer

with open(os.path.join(args.folder, "time_surface.xdmf"), "w") as ofile:
    ofile.write(text)
