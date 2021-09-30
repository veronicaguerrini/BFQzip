#!/usr/bin/env python3

import sys, argparse, os.path, pathlib

Description = """Tool for checking regular FASTQ files

For example
  {exe} INPUT.fastq 
will check for each read, if its DNA and QS sequences 
have same length

--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

########
def checkFASTQ(filename):
    with open(filename, "rt") as F: 
        lines = 0
        q = None
        while q != '': #for each read
            h = F.readline() #header
            d = F.readline() #dna
            p = F.readline() #+
            q = F.readline() #qs
            if(len(d) != len(q)):
                print("{}: {} {}".format(lines*4, len(d), len(q)))
                return False
            lines+=1

    return True

########

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',      help='input file name', type=str, nargs='+')
    args = parser.parse_args()
    if (not check_input(args)):
        return False

    return checkFASTQ(args.input[0])

########

# check correctness of number of input file and define basename for output
def check_input(args):
    ext = args.input[0].split(".")[-1]
    if (ext!="fastq" and ext!="fq"):
        print("Invalid FASTQ!")
        return False

    return True

########

def show_command_line(f):
    f.write("Python command line: ")
    for x in sys.argv:
        f.write(x+" ")
    f.write("\n")

########

if __name__ == '__main__':
    if(main()):
        print("Valid FASTQ file!")
    else:
        print("Invalid FASTQ file!")

