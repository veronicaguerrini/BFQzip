#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, shutil, struct, pathlib

import checkFASTQ as checker

Description = """Tool to compress FASTQ files in internal memory

For example
  {exe} INPUT.fastq -o OUTPUT -1
will produce the files OUTPUT.fq.7z
 
--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

gsufsort_exe = "external/gsufsort/gsufsort"
header_split = "sed -n 1~4p"
qs_split = "sed -n 4~4p"
dna_split = "sed -n 2~4p"
zip7_exe = "7z a -mm=PPMd"
bsc_exe = "external/libbsc/bsc"
spring_reorder_exe = "external/SPRING/build/spring-reorder"

smooth_exe = "src_int_mem/bfq_int"

bwt_ext = ".bwt"
qs_ext  = ".bwt.qs"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',      help='input file name', type=str, nargs='+')
    parser.add_argument('-o','--out', help='output base name (def. input base name)', default="", type=str)  
    parser.add_argument('-T','--mcl', help='minimum context length', default="", type=str)
    parser.add_argument('-Q','--rv',  help='constant replacement value', default="", type=str)
    parser.add_argument('--rebuild',  help='force call step 1', action='store_true',default=False)
    parser.add_argument('--original', help='do not call step 3',action='store_true')
    parser.add_argument('-1', '--m1', help='mode 1: FASTQ', action='store_true',default=True)
    parser.add_argument('-2', '--m2', help='mode 2: DNA+QS', action='store_true')
    parser.add_argument('-3', '--m3', help='mode 3: DNA+QS+H', action='store_true')
    parser.add_argument('-0', '--m0', help='mode 0: do not compress', action='store_true')
    parser.add_argument('--headers',  help='include the headers', action='store_true', default=False)
    parser.add_argument('--reorder',  help='reorder reads (SPRING)', action='store_true', default=False)
    parser.add_argument('-c', '--check', help='Check if the FASTQ is valid', action='store_true', default=False)
    parser.add_argument('-v',         help='verbose: extra info in the log file', default=0, type=int)
    args = parser.parse_args()
    # ---- check number of input files and define basename
    check_input(args)
    define_basename(args)
    # ---- create and open log file
    logfile_name = args.basename + ".log"
    # get main directory
    args.dir = os.path.split(sys.argv[0])[0]
    print("Sending logging messages to file:", logfile_name)

    with open(logfile_name,"w") as logfile:
        
        ##
        if(args.m2): args.m1 = args.m3 = False
        if(args.m3): 
          args.m1 = args.m2 = False
          args.headers = True
        ##

        if(args.m1):
            print(">>> mode 1: FASTQ",file=logfile) 
            print(">>> mode 1: FASTQ") 
        if(args.m2):
            print(">>> mode 2: DNA+QS",file=logfile) 
            print(">>> mode 2: DNA+QS")
        if(args.m3):
            print(">>> mode 3: DNA+QS+H",file=logfile) 
            print(">>> mode 3: DNA+QS+H")

        show_command_line(logfile)
        logfile.flush()

        if(args.reorder):
            if(spring_reorder(args, logfile, logfile_name)==False):
                print("=== ERROR ===")
                print("./spring-reorder not installed (run make SPRING=1)")
                return False

        if len(args.out)==0 : args.out=args.input[0]

        # temporary files
        args.tmp = []
        args.stream = []

        #--- step1: compute BWT+QS
        
        exists = os.path.exists(args.out+bwt_ext) and\
                 os.path.exists(args.out+qs_ext);
           
        if(args.rebuild or not exists):
            if(args.v == 2): print("## STEP 1 ##")
            start = time.time()
            if(step1(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))
        else:
            args.tmp.append(args.out+bwt_ext)
            args.tmp.append(args.out+qs_ext)
        
        exists = os.path.exists(args.out+".h");
        
        if(args.headers and not exists):
            #--- step2: extract headers
            if(args.v == 2): print("## STEP 2 ##")
            start = time.time()
            if(step2(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))
              
        #--- step3: smooth BWT and QS sequences 
        start = time.time()
        if(args.v == 2): print("## STEP 3 ##")
        if(step3(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))

        #--- step4: compute DNA+QS
        if(args.m2 or args.m3):
            if(args.v == 2): print("## STEP 4 ##")
            start = time.time()
            if(step4(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))


        if(not args.m0):#call compressors

            args.output = []
            args.output2 = []
            
            #--- step5: compress
            if(args.v == 2): print("## STEP 5 ##")
            start = time.time()
            print("--- Step 5 ---", file=logfile); logfile.flush()
            if(step5(args, logfile, logfile_name)!=True):
                sys.exit(1)
            if(step5b(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))

            #---- final report
            if(args.v):
                insize = os.path.getsize(args.input[0])

                print("=== results ==="); 
                print("Original:\t{0:.2f} MB".format(insize/(1024*1024)))
                if(args.v == 2):
                    print(args.input[0])
                print("== PPMd ==")
                outsize = 0
                for f in args.output:
                    outsize += os.path.getsize(f)
                print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
                print("Ratio = {0:.2f}".format(outsize/insize))
                if(args.v == 2):
                    for f in args.output:
                       print(f) 
                print("== BSC ==")
                outsize = 0
                for f in args.output2:
                   outsize += os.path.getsize(f)
                print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
                print("Ratio = {0:.2f}".format(outsize/insize))
                if(args.v == 2):
                   for f in args.output2:
                      print(f) 

    return True

##

def step1(args, logfile, logfile_name):
    print("--- Step 1 ---", file=logfile); logfile.flush()
    exe = os.path.join(args.dir, gsufsort_exe)
    options = ""
    if len(args.out)>0 : options+="-o "+args.out
    else : options+=" -o "+args.input[0]
    command = "{exe} {ifile} --bwt --qs {opt}".format(exe=exe, ifile=args.input[0], opt=options)
    print("=== gsufsort ==="); print(command)
    # tmp files
    args.tmp.append(args.basename+".bwt")
    args.tmp.append(args.basename+".bwt.qs")
    return execute_command(command, logfile, logfile_name)

##
def step2(args, logfile, logfile_name):
    print("--- Step 2 ---", file=logfile); logfile.flush()
    ##
    exe = header_split
    ifile = args.input[0]
    ofile = args.out+".h"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== header ===")
    print(command)
    os.system(command)
    if (args.m3): args.stream.append(args.out+".h")
    return True
##

def step3(args, logfile, logfile_name):
    if args.original:
        print("--- Step 3 ---", file=logfile); logfile.flush()
        command = "cp "+ args.input[0] +" "+args.out+".fq" 
        print(command)
        os.system(command)
    else:
        print("--- Step 3 ---", file=logfile); logfile.flush()
        exe = os.path.join(args.dir, smooth_exe)
        options = "-e " + args.tmp[0] + " -q " + args.tmp[1] + " -o "+args.out+".fq"+" -m 5"
        ##additional options
        if len(args.mcl)>0:
            options+=" -k "+args.mcl
        if len(args.rv)>0: 
            options+=" -v "+str(ord(args.rv))
        if(args.headers): #not ignore headers
            options+=" -H "+args.out+".h"
        command = "{exe} {opt}".format(exe=exe, opt=options)
        print("=== smooth-qs ===")
        print(command)
        if(args.m1): args.stream.append(args.out+".fq")
        return execute_command(command, logfile, logfile_name)
    return True
##

def step4(args, logfile, logfile_name):
    print("--- Step 4 ---", file=logfile); logfile.flush()
    exe = dna_split
    ifile = args.out+".fq"
    ofile = args.out+".fq.dna"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== dna ===")
    print(command)
    os.system(command)
    args.stream.append(args.out+".fq.dna")
    ##
    exe = qs_split
    ifile = args.out+".fq"
    ofile = args.out+".fq.qs"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== qs ===")
    print(command)
    os.system(command)
    args.stream.append(args.out+".fq.qs")

    return True

def step5(args, logfile, logfile_name):
    print("--- PPMd ---", file=logfile); logfile.flush()
    exe = zip7_exe
    print("=== PPMd ===")
    for f in args.stream:
        ofile = f+".7z"
        command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=f, ofile=ofile)
        print(command)
        execute_command(command, logfile, logfile_name)
        args.output.append(ofile)
    return True

def step5b(args, logfile, logfile_name):
    print("=== BSC ===", file=logfile); logfile.flush()
    exe = bsc_exe
    print("=== BSC ===")
    for f in args.stream:
        ofile = f+".bsc"
        command = "{exe} e {ifile} {ofile} -T".format(exe=exe, ifile=f, ofile=ofile)
        print(command)
        execute_command(command, logfile, logfile_name)
        args.output2.append(ofile)
    return True

def spring_reorder(args, logfile, logfile_name):
    if(not os.path.exists(spring_reorder_exe)):
        return False
    print("=== SPRING (reorder-only) ===", file=logfile); logfile.flush()
    ##
    exe = spring_reorder_exe
    ifile = args.input[0]
    filename, file_extension = os.path.splitext(ifile)
    ofile = filename+".reordered"+file_extension
    command = "{exe} -i {ifile} -o {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== SPRING (reorder-only) ===")
    print(command)
    execute_command(command, logfile, logfile_name)
    args.input[0] = ofile
    define_basename(args)
    return True

########

# check correctness of number of input file and define basename for output
def check_input(args):
    if args.check:
        print("=== checking FASTQ ==="); 
        if checker.checkFASTQ(args.input[0]):
            print("Valid FASTQ file!")
        else:
            print("Invalid FASTQ file!")
    return True
            
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

# compute hash digest for a file 
def file_digest(name,logfile):
    try:
        hash_command = "{exe} {infile}".format(exe=shasum_exe, infile=name)
        hashsum = subprocess.check_output(hash_command.split(),stderr=logfile)
        hashsum = hashsum.decode("utf-8").split()[0]
    except:
        hashsum = "Error!" 
    return hashsum  

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name):
    try:
        subprocess.check_call(command.split(),stdout=logfile,stderr=logfile)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile_name)
        return False
    return True

def show_command_line(f):
    f.write("Python command line: ") 
    for x in sys.argv:
        f.write(x+" ")
    f.write("\n")   

if __name__ == '__main__':
    main()
