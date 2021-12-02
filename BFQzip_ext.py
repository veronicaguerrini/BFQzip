#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, shutil, struct, pathlib

import checkFASTQ as checker

Description = """Tool to compress FASTQ files in external memory

For example
  {exe} INPUT.fastq -o OUTPUT -1
will produce the files OUTPUT.fq.7z
 
--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

egap_exe = "external/egap/eGap"
header_split = "sed -n 1~4p"
qs_split = "sed -n 4~4p"
dna_split = "sed -n 2~4p"
zip7_exe = "7z a -mm=PPMd"
bsc_exe = "external/libbsc/bsc"

smooth_exe = "src_ext_mem/bfq_ext"

max_read_len = 250 

lcpbytes = "1"
bwt_ext = ".bwt"
qs_ext  = ".bwt.qs"
lcp_ext = "." + lcpbytes +".lcp"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',      help='input file name', type=str, nargs='+')
    parser.add_argument('-o','--out', help='output base name (def. input base name)', default="", type=str)  
    parser.add_argument('-T','--mcl', help='minimum context length', default="", type=str)
    parser.add_argument('-Q','--rv',  help='constant replacement value', default="", type=str)
    parser.add_argument('--rebuild',  help='force call step 1', action='store_true',default=False)
    parser.add_argument('--original', help='do not call step 3', action='store_true')
    parser.add_argument('-1', '--m1', help='mode 1: FASTQ', action='store_true',default=True)
    parser.add_argument('-2', '--m2', help='mode 2: DNA+QS', action='store_true')
    parser.add_argument('-3', '--m3', help='mode 3: DNA+QS+H', action='store_true')
    parser.add_argument('--headers',  help='include the headers', action='store_true', default=False)
    parser.add_argument('-m', '--mem', help='use at most M MBs', default=0, type=int)
    parser.add_argument('-c', '--check', help='Check if the FASTQ is valid', action='store_true', default=True)
    parser.add_argument('-v',         help='verbose: extra info in the log file',action='store_true')
    args = parser.parse_args()
    # ---- check number of input files and define basename
    check_input(args)
    # ---- create and open log file
    logfile_name = args.out + ".log"
    # get main directory
    args.dir = os.path.split(sys.argv[0])[0]
    print("Sending logging messages to file:", logfile_name)

    with open(logfile_name,"w") as logfile:
        
        ##
        if(args.m2): args.m1 = False
        if(args.m3): 
            args.m1 = False
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
        
        if len(args.out)==0 : args.out=args.input[0]

        # temporary files
        args.tmp = []
        args.stream = []

        #--- step1: compute BWT+QS
        
        exists = os.path.exists(args.out+bwt_ext) and\
                 os.path.exists(args.out+qs_ext) and\
                 os.path.exists(args.out+lcp_ext);
                
        if(args.rebuild or not exists):
            start = time.time()
            if(step1(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))
        else:
            args.tmp.append(args.out+bwt_ext)
            args.tmp.append(args.out+qs_ext)
            args.tmp.append(args.out+lcp_ext)
        
        exists = os.path.exists(args.out+".h");
        
        if(args.headers and not exists):
            #--- step2: extract headers
            start = time.time()
            if(step2(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))
    
        #--- step3: smooth BWT and QS sequences 
        start = time.time()
        if(step3(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))

        #--- step4: compute DNA+QS
        if(args.m2 or args.m3):
            start = time.time()
            if(step4(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))

        
        args.output = []
        args.output2 = []

        #--- step5: compress new FASTQ
        start = time.time()
        if(step5(args, logfile, logfile_name)!=True):
            sys.exit(1)
        if(step5b(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))

        # ---- final report
        insize = os.path.getsize(args.input[0])

        print("=== results ==="); 
        print("Original:\t{0:.2f} MB".format(insize/(1024*1024)))
        if(args.v): print(args.input[0])
        print("== PPMd ==")
        outsize = 0
        for f in args.output:
            outsize += os.path.getsize(f)
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))
        if(args.v):
            for f in args.output:
               print(f) 
        print("== BSC ==")
        outsize = 0
        for f in args.output2:
           outsize += os.path.getsize(f)
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))
        if(args.v):
           for f in args.output_b:
              print(f) 

    return True

##

def step1(args, logfile, logfile_name):
    print("--- Step 1 ---", file=logfile); logfile.flush()
    exe = ""
    options = "" 
    command = ""
    if len(args.out)>0 : options+=" -o "+args.out
    else : options+="-o "+args.input[0]
    options += " --lcp --lbytes "+lcpbytes
    exe += os.path.join(args.dir, egap_exe)
    mem_egap = os.path.getsize(args.input[0])//(2**20)
    if(args.mem>0): m = min(mem_egap, args.mem)
    else: m = mem_egap
    command += "{exe} {ifile} --em --mem {mem} --qs {opt}".format(exe=exe, ifile=args.input[0], opt=options, mem=m)
    print("=== egap ==="); print(command)
    # tmp files
    args.tmp.append(args.out+bwt_ext)
    args.tmp.append(args.out+qs_ext)
    args.tmp.append(args.out+lcp_ext)
    return execute_command(command, logfile, logfile_name)

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
    if (args.m3 or args.m4): args.stream.append(args.out+".h")
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
        options = "-e "+args.tmp[0]+" -q "+args.tmp[1]+" -a "+args.tmp[2]+" -o "+args.out+" -l "+str(max_read_len)+" -s 0 -m 5"
        if len(args.mcl)>0:
            options+=" -k "+args.mcl
        if len(args.rv)>0: 
            options+=" -v "+str(ord(args.rv))
        if(args.headers): #include headers
            options+=" -H " + args.out + ".h"
        command = "{exe} {opt}".format(exe=exe, opt=options)
        print("=== smooth-qs ===") #and noise-reduction
        print(command)
        if(args.m1): args.stream.append(args.out+".fq")
        return execute_command(command, logfile, logfile_name)
    return True

##

def step4(args, logfile, logfile_name):
    print("--- Step 4 ---", file=logfile); logfile.flush()
    exe = dna_split
    ifile = args.input[0]
    ofile = args.out+".fq.dna"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== dna ===")
    print(command)
    os.system(command)
    args.stream.append(args.out+".fq.dna")
    ##
    exe = qs_split
    ifile = args.input[0]
    ofile = args.out+".fq.qs"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== qs ===")
    print(command)
    os.system(command)
    args.stream.append(args.out+".fq.qs")

    return True

##

def step5(args, logfile, logfile_name):
    print("--- Step 5 ---", file=logfile); logfile.flush()
    #exe = gzip_exe
    exe = zip7_exe
    print("=== compression ===")
    for f in args.stream:
        #ofile = f+".gz"
        #command = "{exe} {ifile}".format(exe=exe, ifile=f)
        ofile = f+".7z"
        command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=f, ofile=ofile)
        print(command)
        execute_command(command, logfile, logfile_name)
        args.output.append(ofile)
    return True

def step5b(args, logfile, logfile_name):
    print("=== BSC ===", file=logfile); logfile.flush()
    exe = bsc_exe
    for f in args.stream:
        ofile = f+".bsc"
        command = "{exe} e {ifile} {ofile} -T".format(exe=exe, ifile=f, ofile=ofile)
        print(command)
        execute_command(command, logfile, logfile_name)
        args.output2.append(ofile)
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

    if len(args.out)==0:
        args.out=args.input[0]
        args.basename = args.out 
    elif args.out[-1]=="/": 
        pathlib.Path(args.out).mkdir(parents=True, exist_ok=True) 
        tmp = args.input[0].split("/")
        args.basename = args.out+tmp[-1]
        args.out = args.out+tmp[-1]
        args.out = ''.join(args.out.split(".")[0:-1])
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
