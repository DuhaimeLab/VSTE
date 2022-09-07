with open("/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSetToolOutputs/testing_sets_virsorter_output.tsv", 'wt') as outfile:
    outfile.write("index\tcontig\tcategory\n")
    for index in range(1,11):
        filename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSetToolOutputs/Testing_Set_%s/VIRSorter_global-phage-signal.csv" %(index)
        with open(filename, 'rt') as v:
            for line in v:
                c1 = line.split(" ")[0]
                if (c1 == "##"):
                    if ("(" in line):
                        c2 = line.split(" ")[1]
                        category=c2
                if (c1 != "##"):
                    contig=line.split(",")[0]
                    numgenes1=line.split(",")[1]
                    numgenes2=line.split(",")[3]
                    numphagegenes=line.split(",")[5]
                    phagessig=line.split(",")[6]
                    noncaudphagesig=line.split(",")[7]
                    outfile.write('%s\t%s\t%s\n' %(index,contig,category))