#!/usr/bin/env python3

import sys, argparse, os.path, pathlib, random

Description = """Tool for (randomly) reordering reads on FASTQ files

For example
  {exe} INPUT.fastq
gives INPUT.random.fastq with all reads reordered

--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

########

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',      help='input file name', type=str, nargs='+')
    parser.add_argument('-o','--out', help='output base name (def. input base name)', default="", type=str)  
    parser.add_argument('-p', '--paired', help='paired end mode', action='store_true')
    parser.add_argument('-c', '--check', help='Check if the FASTQ is valid', action='store_true', default=False)
    args = parser.parse_args()
    # ---- check number of input files and define basename
    check_input(args, args.input[0])
    if(args.paired):
        if(len(args.input)!=2):
            print("=== ERROR ===")
            print("paired end mode")
            return False
        check_input(args, args.input[1])

    filename, file_extension = os.path.splitext(args.input[0])
    if len(args.out)==0 : 
        if(args.paired):
            filename2, file_extension2 = os.path.splitext(args.input[1])
            ofile = [filename+".random"+file_extension, filename2+".random"+file_extension2]
        else:
            ofile = [filename+".random"+file_extension]
    else:
        if(args.paired):
            ofile = [args.out+"_1"+file_extension, args.out+"_2"+file_extension]
        else:
            ofile = [args.out+file_extension]



    randomized(args, ofile)

########
def randomized(args, ofile):
    
    filename = args.input[0]

    FASTQ = []
    with open(filename, "rt") as F: 
        lines = 0
        q = None
        while q != '': #for each read
            h = F.readline() #header
            d = F.readline() #dna
            p = F.readline() #+
            q = F.readline() #qs
            FASTQ.append((h, d, p, q))
            lines+=1

    if(args.paired):
        filename = args.input[1]
        FASTQ2 = []
        with open(filename, "rt") as F: 
            lines = 0
            q = None
            while q != '': #for each read
                h = F.readline() #header
                d = F.readline() #dna
                p = F.readline() #+
                q = F.readline() #qs
                FASTQ2.append((h, d, p, q))
                lines+=1
    
    #Random sort
    R = [i for i in range(len(FASTQ))]
    random.shuffle(R)

    with open(ofile[0], "w") as F: 
        for t in R:
            F.write(FASTQ[t][0]) #header
            F.write(FASTQ[t][1]) #dna
            F.write(FASTQ[t][2]) #+
            F.write(FASTQ[t][3]) #qs

    if(args.paired):
        with open(ofile[1], "w") as F: 
            for t in R:
                F.write(FASTQ2[t][0]) #header
                F.write(FASTQ2[t][1]) #dna
                F.write(FASTQ2[t][2]) #+
                F.write(FASTQ2[t][3]) #qs


    return True


########

def show_command_line(f):
    f.write("Python command line: ")
    for x in sys.argv:
        f.write(x+" ")
    f.write("\n")

########
# check correctness of number of input file and define basename for output
def check_input(args, ifile):
    if args.check:
        print("=== checking FASTQ ==="); 
        print(ifile)
        if checker.checkFASTQ(ifile):
            print("Valid FASTQ file!")
        else:
            print("Invalid FASTQ file!")
    return True

########

def define_basename(args):
    if len(args.out)==0:
        args.basename = args.input[0]
    elif args.out[-1]=="/": 
        pathlib.Path(args.out).mkdir(parents=True, exist_ok=True) 
        tmp = args.input[0].split("/")
        args.basename = args.out+tmp[-1]
    else:
        args.basename = args.out
    return True

########

if __name__ == '__main__':
    main()

