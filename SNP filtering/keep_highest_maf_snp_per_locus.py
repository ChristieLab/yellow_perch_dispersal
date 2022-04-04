# VCF file
raw_vcf_file = open("raw.maxmiss.7.missindv.7.maf.01.HDplot.recode.vcf", "r")
raw_vcf_array = raw_vcf_file.readlines()
raw_vcf_file.close()



outfile = open("raw.maxmiss.7.missindv.7.maf.01.HDplot.1snp.vcf", "w") # create outfile

# now use array of locus names to loop through locus names in VCF and append matching lines to new file with same headers
#cleaned_vcf_list=[]


good_snp_dict={}
for i in raw_vcf_array:
    # this is trying to recognize the header lines and append to file
    if i.startswith("#"):
        #line = i
        outfile.write(i)
            # at this point we should have header line
    # use the if not # to go to the data
    elif "#" not in i:
        #splits each individual info into respective elements by tabs
        split_ind_line = i.rstrip().split("\t")
        locus=split_ind_line[0]
        base_pair=split_ind_line[1]
        full_line=i
        AF_group=split_ind_line[8]
        split_AF=AF_group.split("=")
        AF=split_AF[2]
        line_for_dict=[AF, full_line]
        if locus not in good_snp_dict.keys():
           good_snp_dict[locus]=line_for_dict
        current_best_AF_line=good_snp_dict[locus]
        if AF > current_best_AF_line[0]:
            good_snp_dict[locus]=line_for_dict

sorted_locus_names=sorted(good_snp_dict.keys(), key=int)
for i in sorted_locus_names:
    dict_entry=good_snp_dict[i]
    outfile.write(dict_entry[1])

outfile.close()