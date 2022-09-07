with open("/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSetToolOutputs/merged.nreuk.kaiju.names.out", 'wt') as outfile:
    outfile.write('Index\tClassified\tContig\tNCBI_taxon\tlen\tID_best\tIDs_all\tSeq\tName\n')
    for index in range(1,11):
        filename="/scratch/duhaimem_root/duhaimem/shared_data/VSTE/TestingSetToolOutputs/Testing_Set_%s/%s.nreuk.kaiju.names.out" %(index,index)
        with open(filename, 'rt') as v:
            for line in v:
                outfile.write('%s\t%s' %(index, line))