#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, shutil, struct, pathlib, threading, logging

import checkFASTQ as checker

Description = """Tool to compress FASTQ files in internal memory

For example
  {exe} INPUT.fastq -o OUTPUT
will produce the files OUTPUT.fq.7z and OUTPUT.fq.bsc
 
--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

bfqzip_exe = "python3 BFQzip.py"
zip7_exe = "7z a -mm=PPMd"
bsc_exe = "external/libbsc/bsc"
spring_reorder_exe = "external/SPRING/build/spring-reorder"
random_reorder_py= "randomFASTQ.py"
random_reorder_exe = "python3 randomFASTQ.py"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',      help='input file name', type=str, nargs='+')
    parser.add_argument('-o','--out', help='output base name (def. input base name)', default="", type=str)  
    parser.add_argument('-T','--mcl', help='minimum context length', default="", type=str)
    parser.add_argument('-Q','--rv',  help='constant replacement value', default="", type=str)
    parser.add_argument('-H', '--headers', help='store original headers', action='store_true')
    parser.add_argument('-0', '--m0', help='mode 0: do not compress', action='store_true')
    parser.add_argument('-p', '--paired', help='paired end mode', action='store_true')
    parser.add_argument('-t','--threads',  help='multithreading', default=0, type=int)
    parser.add_argument('--reorder',  help='reorder reads (0: no reorder, 1: random, 2: SPRING)', default=0, type=int)
    parser.add_argument('-c', '--check', help='Check if the FASTQ is valid', action='store_true', default=False)
    parser.add_argument('-v',         help='verbose: extra info in the log file', default=0, type=int)
    args = parser.parse_args()
    # ---- check number of input files and define basename
    check_input(args, args.input[0])
    if(args.paired):
        if(len(args.input)!=2):
            print("=== ERROR ===")
            print("paired end mode")
            return False
        check_input(args, args.input[1])
    define_basename(args)
    # ---- create and open log file
    logfile_name = args.basename + ".log"
    print("Sending logging messages to file:", logfile_name)

    print("Threads =", args.threads)
    with open(logfile_name,"w") as logfile:

        ##number of threads
        show_command_line(logfile)
        logfile.flush()

        if(args.reorder>0):
            start = time.time()
            reorder = True
            #random reorder
            if(args.reorder==1):
                if(random_reorder(args, logfile, logfile_name)==False):
                    return False
            #SPRING-reorder
            elif(args.reorder==2):
                if(spring_reorder(args, args.input[0], logfile, logfile_name)==False):
                    reorder = False
                if(args.paired):
                    if(spring_reorder(args, args.input[1], logfile, logfile_name)==False):
                        reorder = False
                if(not reorder):
                    return False
            print("Elapsed time: {0:.4f}".format(time.time()-start))

        start = time.time()
        blocks = []
        split_fastq(args, args.input[0], blocks, logfile, logfile_name)
        if(args.paired):
            split_fastq_2(args, args.input[0], args.input[1], logfile, logfile_name)
        print("Elapsed time: {0:.4f}".format(time.time()-start))

        files=""
        for i in blocks:
            filename = i[0]
            files += filename+".fq "

        ####
        threads = list()

        options = "--rebuild -0 -v 0"
        if args.headers:
            options+=" --headers"
        if args.check:
            options+=" --check"
        if len(args.mcl)>0:
            options+=" -T "+args.mcl
        if len(args.rv)>0: 
            options+=" -Q "+args.rv

        start = time.time()
        if(args.threads==0):
            for i in blocks:
                filename = i[0]
                bfqzip(filename, options, logfile, logfile_name)

        #multithreading
        else:
            for i in blocks:
                filename = i[0]
                x = threading.Thread(target=bfqzip, args=(filename, options, logfile, logfile_name,))
                threads.append(x)
                x.start()

            for index, thread in enumerate(threads):
                logging.info("Main    : before joining thread %d.", index)
                thread.join()
                logging.info("Main    : thread %d done", index)

        #remove files
        for i in blocks:
            filename = i[0]
            os.remove(filename)
            os.remove(filename+".bwt")
            os.remove(filename+".bwt.qs")
            if args.headers:
                os.remove(filename+".h")
            os.remove(filename+".log")
            print(i[0], i[1], file=logfile)

        print("Total time: {0:.4f}".format(time.time()-start))

        print(args.out)

        print()
        print("=== MERGE ===")
        #concatenate files
        filename, file_extension = os.path.splitext(args.input[0])

        if len(args.out)==0: 
            if(args.paired):
                filename2, file_extension2 = os.path.splitext(args.input[1])
                ofile = [filename+".cat"+file_extension, filename2+".cat"+file_extension2]
            else:
                ofile = [filename+".cat"+file_extension]
        else:
            if(args.paired):
                ofile = [args.out+"_1"+file_extension, args.out+"_2"+file_extension]
            else:
                ofile = [args.out+file_extension]
        
        if(args.paired):
            if(os.path.exists(ofile[0])):
                os.remove(ofile[0])
            if(os.path.exists(ofile[1])):
                os.remove(ofile[1])
            for i in blocks:
                ifile=i[0]+".fq"
                with open(ifile, 'r') as f_in:
                    num_lines_1 = i[1]
                    
                    with open(ofile[0], 'a') as f_out:
                        num_lines = 0
                        for line in f_in:
                            f_out.write(line)
                            num_lines+=1
                            if(num_lines==num_lines_1):
                                break
                    with open(ofile[1], 'a') as f_out:
                        for line in f_in:
                            f_out.write(line)
        else:
            cmd = "cat "+files+"> "+ofile[0]
            print(cmd)
            start = time.time()
            os.system(cmd)

        print("Elapsed time: {0:.4f}".format(time.time()-start))

        #remove files
        for i in blocks:
            filename = i[0]
            os.remove(filename+".fq")

        args.stream = []
        args.stream.extend(ofile)
        #---- final report

        if(not args.m0):#call compressors

            args.output = []
            args.output2 = []
            
            #--- step5: compress
            if(args.v == 2): print("## STEP 5 ##")
            print("--- Step 5 ---", file=logfile); logfile.flush()
            print()
            start = time.time()
            if(args.threads==0):
                for f in args.stream:
                    if(step5(args, f, logfile, logfile_name)!=True):
                        sys.exit(1)
            else:
                threads = list()
                for f in args.stream:
                    x = threading.Thread(target=step5, args=(args, f, logfile, logfile_name,))
                    threads.append(x)
                    x.start()
                for index, thread in enumerate(threads):
                    logging.info("Main    : before joining thread %d.", index)
                    thread.join()
                    logging.info("Main    : thread %d done", index)



            print("Elapsed time: {0:.4f}".format(time.time()-start))
            start = time.time()

            if(args.threads==0):
                for f in args.stream:
                    if(step5b(args, f, logfile, logfile_name)!=True):
                        sys.exit(1)
            else:
                threads = list()
                for f in args.stream:
                    x = threading.Thread(target=step5b, args=(args, f, logfile, logfile_name,))
                    threads.append(x)
                    x.start()
                for index, thread in enumerate(threads):
                    logging.info("Main    : before joining thread %d.", index)
                    thread.join()
                    logging.info("Main    : thread %d done", index)


            print("Elapsed time: {0:.4f}".format(time.time()-start))

            #---- final report
            if(args.v):
                insize = os.path.getsize(args.input[0])

                print()
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
def thread_function(name):
    logging.info("Thread %s: starting", name)
    time.sleep(3)
    logging.info("Thread %s: finishing", name)

##

def bfqzip(ifile, options, logfile, logfile_name):
    print("--- BFQzip ---", file=logfile); logfile.flush()
    ##
    exe = bfqzip_exe
    command = "{exe} {ifile} {options}".format(exe=exe, ifile=ifile, options=options)
    print("\n=== BFQzip ===\n")
    print(command)
    os.system(command)
    return True
##

def split_fastq(args, ifile, blocks, logfile, logfile_name):
    print("=== Split FASTQ ===", file=logfile)
    print("=== Split FASTQ ===")
    ##
    filename, file_extension = os.path.splitext(ifile)


    with open(ifile, 'r') as f_in:
        num_lines = sum(1 for line in f_in)
        f_in.seek(0)

        num_reads = num_lines//4

        if(args.threads==0):
            size_block = num_reads
        else:
            size_block = num_reads//args.threads
        num_blocks = num_reads//size_block

        for i in range(num_blocks):
            ofile = filename+"."+str(i+1)+file_extension
            with open(ofile, 'w') as f_out:
                num_lines = 0
                num_reads = 0
                for line in f_in:
                    f_out.write(line)
                    num_lines+=1
                    if(num_lines%4==0): 
                        num_reads+=1
                    if(num_reads==size_block and i+1 != num_blocks):
                        break
            blocks.append((ofile, num_lines))

    print('{} blocks of ~{} lines:'.format(num_blocks,size_block))

    return True

def split_fastq_2(args, ifile1, ifile2, logfile, logfile_name):
    print("=== Split FASTQ 2 ===", file=logfile)
    print("=== Split FASTQ 2 ===")
    ##
    filename, file_extension = os.path.splitext(ifile1)


    with open(ifile2, 'r') as f_in:
        num_lines = sum(1 for line in f_in)
        f_in.seek(0)

        num_reads = num_lines//4

        if(args.threads==0):
            size_block = num_reads
        else:
            size_block = num_reads//args.threads
        num_blocks = num_reads//size_block

        for i in range(num_blocks):
            ofile = filename+"."+str(i+1)+file_extension
            with open(ofile, 'a') as f_out:
                num_lines = 0
                num_reads = 0
                for line in f_in:
                    f_out.write(line)
                    num_lines+=1
                    if(num_lines%4==0): 
                        num_reads+=1
                    if(num_reads==size_block and i+1 != num_blocks):
                        break
            #blocks.append((ofile, num_lines))

    print('{} blocks of ~{} lines:'.format(num_blocks,size_block))

    return True

########

def step5(args, f, logfile, logfile_name):
    print("--- PPMd ---", file=logfile); logfile.flush()
    exe = zip7_exe
    print("=== PPMd ===")
    ofile = f+".7z"
    command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=f, ofile=ofile)
    print(command)
    execute_command(command, logfile, logfile_name)
    args.output.append(ofile)
    return True

def step5b(args, f, logfile, logfile_name):
    print("=== BSC ===", file=logfile); logfile.flush()
    exe = bsc_exe
    print("=== BSC ===")
    ofile = f+".bsc"
    command = "{exe} e {ifile} {ofile} -b64".format(exe=exe, ifile=f, ofile=ofile)
    print(command)
    execute_command(command, logfile, logfile_name)
    args.output2.append(ofile)
    return True


########

def spring_reorder(args, ifile, logfile, logfile_name):
    if(not os.path.exists(spring_reorder_exe)):
        print("=== ERROR ===")
        print("./spring-reorder not installed (run make SPRING=1)")
        return False
    print("=== SPRING (reorder-only) ===", file=logfile); logfile.flush()
    ##
    exe = spring_reorder_exe
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

def random_reorder(args, logfile, logfile_name):
    ifile = args.input[0]
    if(not os.path.exists(random_reorder_py)):
        print("=== ERROR ===")
        print("{} not found".format(random_reorder_py))
        return False
    print("=== RANDOM (reorder) ===", file=logfile); logfile.flush()
    ##
    exe = random_reorder_exe
    filename, file_extension = os.path.splitext(ifile)
    if(args.paired):
        filename2, file_extension2 = os.path.splitext(args.input[1])
        ofile = [filename+".random"+file_extension, filename2+".random"+file_extension2]
    else:
        ofile = [filename+".random"+file_extension]
    ##
    if(args.paired):
        command = "{exe} {ifile1} {ifile2} --paired".format(exe=exe, ifile1=args.input[0], ifile2=args.input[1])
    else:
        command = "{exe} {ifile}".format(exe=exe, ifile=ifile)
    print("=== RANDOM (reorder) ===")
    print(command)
    execute_command(command, logfile, logfile_name)
    args.input[0] = ofile[0]
    if(args.paired):
        args.input[1] = ofile[1]
    define_basename(args)
    return True

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
