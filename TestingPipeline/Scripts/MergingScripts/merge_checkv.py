with open("/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSetToolOutputs/testing_sets_checkv_output.tsv", 'wt') as outfile:
    j=0
    for index in range(1,11):
        filename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSetToolOutputs/Testing_Set_%s/quality_summary.tsv" %(index)
        with open(filename, 'rt') as v:
            if (j==0):
                i=0
                for line in v:
                    if (i==0):
                        outfile.write('Index\t%s' %(line))
                    if (i>0):
                        outfile.write('%s\t%s' %(index, line))
                    i+=1
            if (j>0):
                i=0
                for line in v:
                    if (i>0):
                        outfile.write('%s\t%s' %(index, line))
                    i+=1
        j+=1