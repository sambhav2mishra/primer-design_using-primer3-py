import primer3

seq_dict = {}
no_primer = []
one_set = []
two_set = []
params_dict = {
    'PRIMER_OPT_SIZE': 20,  
    'PRIMER_MIN_SIZE': 18,  
    'PRIMER_MAX_SIZE': 27,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 55.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_PRODUCT_SIZE_RANGE': [[120,200]]
}


with open('sequences.fa','r') as seq:
    blocks = seq.read().split('>')
    for block in blocks:
        if block:
            segments = block.splitlines()
            #print(segments)
            seq_dict['SEQUENCE_ID'] =  segments[0]
            seq_dict['SEQUENCE_TEMPLATE'] = ''.join(segments[1:len(segments)])
            primers = primer3.bindings.designPrimers(seq_dict,params_dict)
            #print(seq_dict)
           
            if ("PRIMER_LEFT_0_SEQUENCE" in primers) and ("PRIMER_LEFT_1_SEQUENCE" in primers):
                data = segments[0]+"\t"+primers["PRIMER_LEFT_0_SEQUENCE"]+","+primers["PRIMER_RIGHT_0_SEQUENCE"]+"\t"+str(primers["PRIMER_PAIR_0_PRODUCT_SIZE"])+"\t"+str(primers["PRIMER_LEFT_0_TM"])+"\t"+str(primers["PRIMER_RIGHT_0_TM"])+"\t"+primers["PRIMER_LEFT_1_SEQUENCE"]+","+primers["PRIMER_RIGHT_1_SEQUENCE"]+"\t"+str(primers["PRIMER_PAIR_1_PRODUCT_SIZE"])+"\t"+str(primers["PRIMER_LEFT_1_TM"])+"\t"+str(primers["PRIMER_RIGHT_1_TM"])+"\n"
                two_set.append(data)
            elif ("PRIMER_LEFT_0_SEQUENCE" in primers):
                data1 = segments[0]+"\t"+primers["PRIMER_LEFT_0_SEQUENCE"]+","+primers["PRIMER_RIGHT_0_SEQUENCE"]+"\t"+str(primers["PRIMER_PAIR_0_PRODUCT_SIZE"])+"\t"+str(primers["PRIMER_LEFT_0_TM"])+"\t"+str(primers["PRIMER_RIGHT_0_TM"])+"\n"
                one_set.append(data1)
            else:
                no_primer.append(segments[0])

with open("one_set_seq.txt", 'w') as fo:
    for line in one_set:
        fo.write(line)
    
with open("two_set_seq.txt", 'w') as fo:
    for line in two_set:
        fo.write(line)
        
with open("No_primer.txt", "w") as fo:
    for line in no_primer:
        fo.write(line)

