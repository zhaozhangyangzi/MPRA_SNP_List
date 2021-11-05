import sys
import re

if len(sys.argv)<2:
	print 'Usage: python process_LD_filter.py input_file output_file'
	sys.exit()
else:
	infileName=sys.argv[1]
	outfileName=sys.argv[2]

# read in the files
infile=open(infileName,'rU') #open the input file
LineNumber1=0;
SnpData={}
for temp in infile:
	SnpData[LineNumber1]=temp.strip('\n')
	LineNumber1+=1
infile.close()

outfile=open(outfileName,'w') 
LineNumber2=0
for i in range(0,len(SnpData)):
	Line=SnpData[i].strip('\n')
	Line=Line.split('\t')
	lead_Snps=Line[0]
	DL_Snps=Line[1]
	cohort=Line[2]
	DL_Snps=DL_Snps.split(';')
	DL_value=[float(x.split(',')[1]) for x in DL_Snps]
	SNP_name=[x.split(',')[0] for x in DL_Snps]
	DL_SNP_dict=dict(zip(SNP_name, DL_value))
	DL_Snps_sig=[i for i,j in DL_SNP_dict.items() if j >= 0.8]
	DL_r2_sig=[j for j,j in DL_SNP_dict.items() if j >= 0.8]
	if len(DL_r2_sig)==0:
		DL_Snps_sig1=lead_Snps
		DL_r2_sig1='NA'
		causalSNP=DL_Snps_sig1+'\t'+DL_r2_sig1+'\t'+lead_Snps+'\t'+cohort+'\n'
		outfile.write(causalSNP+'\n')
	else:
		for m in range(0,len(DL_Snps_sig)):
			causalSNP=DL_Snps_sig[m]+'\t'+str(DL_r2_sig[m])+'\t'+lead_Snps+'\t'+cohort+'\n'
			outfile.write(causalSNP+'\n')
	LineNumber2+=1
outfile.close()


