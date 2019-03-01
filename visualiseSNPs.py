__author__ = 'Samantha Campbell'

import argparse
import vcf
import re


def getSampleNames(vcf_reader):

    sample_info = vcf_reader.samples
    output_names = []

    for sample in sample_info:
        sample_name = re.split(r'[=\s]\s*', sample)
        output_names.append(sample_name[0])
    return output_names


def mergeReplicates(samples):

    conditions = {'0': [], '5': [], '10': [], '16': [], '20': [], '29': [], '1': [], 'control': []}

    for sample in samples:
        # sort samples by SP no.
        if 'SP' in sample:
            sp = re.split(r'[-\s]\s*', sample)
            num = re.sub("\D", "", sp[1])
        elif 'sp' in sample:
            sp = re.split(r'[_\s]\s*', sample)
            num = re.split(r'[sp\s]\s*', sp[0])
            num = num[3]
        else:
            num = "control"
        conditions.setdefault(num, []).append(sample)

    return conditions



def writeSampleFiles(vcf_reader, samples):
    file_handle = {}
    for sample in samples:
        outfile = sample + ".bed"
        file_handle[sample] = open(outfile, 'w')
        print ("Opened handle for {}".format(sample))


    # for each sample, write each entry to a new file with the data from that column


    for record in vcf_reader:
        for sample in samples:
            print ("It gets here...")
            #   filter on the variant for each sample - only include in the output if variant at that position
            genotype = str(record.genotype(sample)['GT'])
            print (genotype)
            if genotype != "0/0":
                chr = str(record.CHROM)
                start = str(record.start)
                end = str(record.end)
                print ("added:  " + genotype)
                file_handle[sample].write(chr + "\t" + start + "\t" + end + "\n")

    for sample in samples:
        file_handle[sample].close()
        print("Closed file for {}".format(sample))


def writeMergedSampleFiles(vcf_reader, conditions):

    merged = []
    for condition in conditions:
        outfile = "Lmex_SP" + condition + ".bed"
        merged.append(outfile)
        file_out = open(outfile, 'w+')
        originals = conditions[condition]
       # print (originals)

        for file in originals:
            #print (file)
            for record in vcf_reader:
                genotype = str(record.genotype(file)['GT'])
                if genotype != "0/0":
                    chr = str(record.CHROM)
                    start = str(record.start)
                    end = str(record.end)
                    #print(genotype)

                    file_out.write(chr + "\t" + start + "\t" + end + "\n")
        file_out.close()
    return merged


def writeBasicConfig(configName):

    config = open(configName, 'w+')
    config.write(
        "# circos.conf \n\nkaryotype = data/karyotype/karyotype.lmexicana.txt \nchromosomes_display_default = yes\n\n")
    config.write("<fonts> \n<<include etc/fonts.conf>> \n</fonts>\n\n")
    config.write("<colors> \n<<include etc/colors.conf>> \n</colors> \n\n")
    config.write("<ideogram>\n\n<spacing>\ndefault = 0.004r\n</spacing>\n\nradius    = 0.3r\nthickness = 20p\n"
                 "fill      = yes\nshow_label       = no\nlabel_font       = default\nlabel_radius     = 1r + 75p\n"
                 "label_size       = 30\nlabel_parallel   = yes\n\n</ideogram>\n\n")
    config.close()

def completeConfig(configName):

    required = "################ \n\n<image>\n<<include etc/image.conf>>\n</image>\n\n" \
               "<<include etc/colors_fonts_patterns.conf>>\n\n<<include etc/housekeeping.conf>>"
    with open(configName, "a") as myfile:
        myfile.write(required)



def writePlots(configName, samples):

    tracks = len(samples)

    colours = ["blue", "dblue", "lblue", "red", "green", "purple", "dred", "orange", "dorange", "lred", "dgreen", "lgreen"]

    if tracks < len(colours):

        with open(configName, "a") as myfile:
            myfile.write("<highlights>\n\n")

            myfile.write("<highlight>\n\nfile\t= /Users/samanthacampbell/PycharmProjects/vcfToCircos/Lmex_SP16.bed \nstroke_color = blue\nfill_color = dblue\nr0\t= 1.0r\nr1\t= 2.0r\n\n</highlight>\n\n")

            myfile.write("</highlights>\n\n")



#################   main code   ##################

parser = argparse.ArgumentParser(description='Read vcf file')
parser.add_argument('--vcfIn', required=True, type=str, help='vcf input file')
args = parser.parse_args()

#   open vcf file
file = args.vcfIn
vcf_reader = vcf.Reader(open(file, 'r'))


#   get sample names
output_samples = getSampleNames(vcf_reader)

#   merge replicates
merged_samples = mergeReplicates(output_samples)

#   write variant locations to a bed file for each sample
#writeSampleFiles(vcf_reader, output_samples)
#merged_files = writeMergedSampleFiles(vcf_reader, merged_samples)   ######  THIS ISN'T WORKING!!

print ("### Variants written to file    ###")

#   write basic config file
configFile = "/Users/samanthacampbell/Downloads/circos-0.69-2/circos.lmexicana.conf"
#writeBasicConfig(configFile)

print ("### Writing config file...  ###")

#   write plots to file
#writePlots(configFile, merged_files)

completeConfig(configFile)
