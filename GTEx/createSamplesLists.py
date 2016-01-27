#### This script converts the training and testing sample files into a different format that can be manipulated easily in R for fitting the regression models
#### Create dictionaries with training and testing sample IDs ####
train_file = open("/ifs/scratch/c2b2/ip_lab/ak3808/GTEx_processed/phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train" , 'r')
train_dict = train_file.readlines()
test_file = open("/ifs/scratch/c2b2/ip_lab/ak3808/GTEx_processed/phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test",'r')
test_dict = test_file.readlines()

train_dict = map(lambda x: x.split() , train_dict)
test_dict = map(lambda x: x.split(), test_dict)

temp = {}
for x in train_dict:
	temp[x[0]] = x[1:]
train_dict = temp

temp = {}
for x in test_dict:
	temp[x[0]] =  x[1:]
test_dict = temp	
temp = {}


####

f1 = open("../tmp/train_samples.txt",'w')
for k in train_dict.keys():
	for v in train_dict[k]:
		f1.write(str(v)+"\t"+str(k)+"\n")


f2 = open("../tmp/test_samples.txt",'w')
for k in test_dict.keys():
        for v in test_dict[k]:
                f2.write(str(v)+"\t"+str(k)+"\n")
