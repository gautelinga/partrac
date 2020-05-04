import argparse
from lxml import etree
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert xdmf to list of times and files")
    parser.add_argument("filename", type=str, help="XDMF file")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    try:
        tree = etree.parse(args.filename)
    except:
        name = os.path.basename(args.filename).split(".")
        newfilename = os.path.join(
            os.path.dirname(args.filename),
            ".".join(name[:-1]) + "_mod." + name[-1])
        # with open(args.filename)
        print("Making new file: " + newfilename)
        with open(args.filename, "r") as fin:
            content = fin.read()
        content += "    </Grid>\n</Domain>\n</Xdmf>"
        with open(newfilename, "w") as fout:
            fout.write(content)
        tree = etree.parse(newfilename)

    root = tree.getroot()

    it = 0
    entries = []
    for el in root[0]:
        if el.tag == "Grid":
            for subel in el:
                if subel.tag == "Time":
                    times = [int(timestr) for timestr in
                             subel[0].text.strip().split(" ")]
                if subel.tag == "Grid":
                    entries_loc = dict()
                    for subsubel in subel:
                        if subsubel.tag == "Attribute":
                            entries_loc[subsubel.attrib["Name"]] = \
                                subsubel[0].text.strip()
                    entries.append(entries_loc)
                    it = it + 1

    dirname = os.path.dirname(args.filename)
    ofile = open(os.path.join(dirname, "timestamps.dat"), "w")
    req_fields = ["density", "pressure", "u_x", "u_y", "u_z"]
    for t, entries_loc in zip(times, entries):
        for req_field in req_fields:
            assert(req_field in entries_loc.keys())
        files = [entries_loc[req_field].split(":")[0] for req_field in req_fields]
        file = files[0]
        assert(all([files[i] == file for i in range(1, len(files))]))
        assert(entries_loc["is_solid"].split(":")[0] == "output_is_solid.h5")
        # print(t, file)
        ofile.write("{}\t{}\n".format(t, file))
    ofile.close()


if __name__ == "__main__":
    main()

    


