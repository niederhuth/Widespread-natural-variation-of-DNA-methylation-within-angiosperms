import sys
import pandas
import pybedtools
import math
import numpy

#converts allc file to special bed file format for use with pybedtools
#use: allc2bed("allc file")
def allc2bed(allc):
    a = pandas.read_table(allc)
    a['pos2'] = a.pos
    a['name'] = a.index 
    a['score'] = "."
    a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
    a = pybedtools.BedTool.from_dataframe(a)
    return a

#map cytosines from allc to features in bed file
#use: map_allc2feature( "allc file", "bed file of features")
#outputs tab-delimited file: column 1 = feature name, column 2 = sequence context, column 3 = methylated reads, 
#column 4 = total reads, column 5 = is site methylated or unmethylated (determined by methylpy binomial test, 1 = yes, 0 = no)
def map_allc2feature(allc,features_bed):
    mC_bed = allc2bed(allc)
    f_bed = pybedtools.BedTool(features_bed)
    mapping = pybedtools.bedtool.BedTool.intersect(mC_bed,f_bed,wa=True,wb=True)
    m = pandas.read_table(mapping.fn, header=None, usecols = [13,6,7,8,9])
    m = m[[13,6,7,8,9]]
    m = m.sort_values(by = 13)
    return m

#calculate the genome-wide weighted methylation using plant specific contexts
#use: weighted_mC_plant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = weighted methylation (decimal), column 5 = weighted methylation (as %)
def weighted_mC_plant(allc, mC_results, cutoff=0):
    CG = mCG = CHG = mCHG = CHH = mCHH = CNN = mCNN = 0
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + int(c[5])
                    mCG = mCG + int(c[4])
                elif c[3].endswith("G"):
                    CHG = CHG + int(c[5])
                    mCHG = mCHG + int(c[4])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + int(c[5])
                    mCNN = mCNN + int(c[4])
                else:
                    CHH = CHH + int(c[5])
                    mCHH = mCHH + int(c[4])     
    with open(mC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Weighted_mC","%Weighted_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CHG",str(CHG),str(mCHG),str(numpy.float64(mCHG)/numpy.float64(CHG)),str(pCHG*100)]) + '\n')
        out.write('\t'.join(["CHH",str(CHH),str(mCHH),str(numpy.float64(mCHH)/numpy.float64(CHH)),str(pCHH*100)]) + '\n')
        out.write('\t'.join(["CNN",str(CNN),str(mCNN),str(numpy.float64(mCNN)/numpy.float64(CNN)),str(pCNN*100)]) + '\n')

#calculate the genome-wide weighted methylation using non-plant specific contexts
#use: weighted_mC_nonplant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = weighted methylation (decimal), column 5 = weighted methylation (as %)
def weighted_mC_nonplant(allc, mC_results, cutoff=0):    
    CG = mCG = CH = mCH = CNN = mCNN = 0    
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + int(c[5])
                    mCG = mCG + int(c[4])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + int(c[5])
                    mCNN = mCNN + int(c[4])
                else:
                    CH = CH + int(c[5])
                    mCH = mCH + int(c[4])      
    with open(mC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Weighted_mC","%Weighted_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CH",str(CH),str(mCH),str(numpy.float64(mCH)/numpy.float64(CH)),str(pCH*100)]) + '\n')
        out.write('\t'.join(["CN",str(CNN),str(mCNN),str(numpy.float64(mCNN)/numpy.float64(CNN)),str(pCNN*100)]) + '\n')

#calculate the genome-wide percent methylated sites using plant specific contexts. 
#This relies on the binomial test conducted by methylpy to determine if a site is methylated or unmethylated
#use: weighted_mC_plant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = percent methylation (decimal), column 5 = percent methylation (as %)
def percent_mC_plant(allc, pC_results, cutoff=0):
    CG = mCG = CHG = mCHG = CHH = mCHH = CNN = mCNN = 0
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + 1
                    mCG = mCG + int(c[6])
                elif c[3].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[6])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + 1
                    mCNN = mCNN + int(c[6])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[6])       
    with open(pC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Percent_mC","%Percent_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CHG",str(CHG),str(mCHG),str(numpy.float64(mCHG)/numpy.float64(CHG)),str(pCHG*100)]) + '\n')
        out.write('\t'.join(["CHH",str(CHH),str(mCHH),str(numpy.float64(mCHH)/numpy.float64(CHH)),str(pCHH*100)]) + '\n')
        out.write('\t'.join(["CNN",str(CNN),str(mCNN),str(numpy.float64(mCNN)/numpy.float64(CNN)),str(pCNN*100)]) + '\n')

#calculate the genome-wide percent methylated sites using nonplant specific contexts. 
#This relies on the binomial test conducted by methylpy to determine if a site is methylated or unmethylated
#use: weighted_mC_plant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = percent methylation (decimal), column 5 = percent methylation (as %)        
def percent_mC_nonplant(allc, pC_results, cutoff=0):
    CG = mCG = CH = mCH = CNN = mCNN = 0    
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + 1
                    mCG = mCG + int(c[6])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + 1
                    mCNN = mCNN + int(c[6])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[6])         
    with open(mC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Percent_mC","%Percent_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CH",str(CH),str(mCH),str(numpy.float64(mCH)/numpy.float64(CH)),str(pCH*100)]) + '\n')
        out.write('\t'.join(["CN",str(CNN),str(mCNN),str(numpy.float64(mCNN)/numpy.float64(CNN)),str(pCNN*100)]) + '\n')
        
#calculate methylation levels for a region using plant specific contexts
#use: calc_feature_mC_plant( <input file>, <output file>, <cutoff for coverage needed to include site, default=0>)
#takes as input, output from map_allc2feature()
#outputs tab-delimited file: column 1 = feature name, column 2 = total CG sites, column 3 = methylated CG sites, 
#column 4 = percent methylated CGs, column 5 = total CG reads, column 6 = methylated CG reads, 
#column 7 = weighted CG methylation, column 8 = total CHG sites, column 9 = methylated CHG sites, 
#column 10 = percent methylated CHGs, column 11 = total CHG reads, column 12 = methylated CHG reads, 
#column 13 = weighted CHG methylation, column 14 = total CHH sites, column 15 = methylated CHH sites, 
#column 16 = percent methylated CHHs, column 17 = total CHH reads, column 18 = methylated CHH reads, 
#column 19 = weighted CHH methylation 
def calc_feature_mC_plant( allc2feature, feature_mC, cutoff=0):
    CG = mCG = CGr = mCGr = CHG = mCHG = CHGr = mCHGr = CHH = mCHH = CHHr = mCHHr = 0
    name = "none"
    out = open(feature_mC, "w")
    for c in allc2feature.itertuples():
        if name == "none":
            name = c[1]
            out.write('\t'.join(["Feature","Num_CG","mCG","Percent_mCG","Num_CG_reads","mCG_reads","Weighted_mCG","Num_CHG","mCHG","Percent_mCHG","Num_CHG_reads","mCHG_reads","Weighted_mCHG","Num_CHH","mCHH","Percent_mCHH","Num_CHH_reads","mCHH_reads","Weighted_mCHH"]) + '\n')
        elif c[1] != name:
            out.write('\t'.join([str(name),str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(CGr),str(mCGr),str(numpy.float64(mCGr)/numpy.float64(CGr)),str(CHG),str(mCHG),str(numpy.float64(mCHG)/numpy.float64(CHG)),str(CHGr),str(mCHGr),str(numpy.float64(mCHGr)/numpy.float64(CHGr)),str(CHH),str(mCHH),str(numpy.float64(mCHH)/numpy.float64(CHH)),str(CHHr),str(mCHHr),str(numpy.float64(mCHHr)/numpy.float64(CHHr))]) + '\n')
            name = c[1]
            CG = mCG = CGr = mCGr = CHG = mCHG = CHGr = mCHGr = CHH = mCHH = CHHr = mCHHr = 0
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                elif c[2].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[5])
                    CHGr = CHGr + int(c[4])
                    mCHGr = mCHGr + int(c[3])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[5])
                    CHHr = CHHr + int(c[4])
                    mCHHr = mCHHr + int(c[3])
        elif c[1] == name:
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                elif c[2].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[5])
                    CHGr = CHGr + int(c[4])
                    mCHGr = mCHGr + int(c[3])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[5])
                    CHHr = CHHr + int(c[4])
                    mCHHr = mCHHr + int(c[3])
    out.write('\t'.join([str(name),str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(CGr),str(mCGr),str(numpy.float64(mCGr)/numpy.float64(CGr)),str(CHG),str(mCHG),str(numpy.float64(mCHG)/numpy.float64(CHG)),str(CHGr),str(mCHGr),str(numpy.float64(mCHGr)/numpy.float64(CHGr)),str(CHH),str(mCHH),str(numpy.float64(mCHH)/numpy.float64(CHH)),str(CHHr),str(mCHHr),str(numpy.float64(mCHHr)/numpy.float64(CHHr))]) + '\n')
    out.close()
              
#calculate methylation levels for a region using nonplant specific contexts
#use: calc_feature_mC_nonplant( <input file>, <output file>, <cutoff for coverage needed to include site, default=0>)
#takes as input, output from map_allc2feature()
#outputs tab-delimited file: column 1 = feature name, column 2 = total CG sites, column 3 = methylated CG sites, 
#column 4 = percent methylated CGs, column 5 = total CG reads, column 6 = methylated CG reads, 
#column 7 = weighted CG methylation, column 8 = total CH sites, column 9 = methylated CH sites, 
#column 10 = percent methylated CHs, column 11 = total CH reads, column 12 = methylated CH reads, 
#column 13 = weighted CH methylation
def calc_feature_mC_nonplant( allc2feature, feature_mC, cutoff=0):
    CG = mCG = CGr = mCGr = CH = mCH = CHr = mCHr = 0
    name = "none"
    out = open(feature_mC, "w")
    for c in allc2feature.itertuples():
        if name == "none":
            name = c[1]
            out.write('\t'.join(["Feature","Num_CG","mCG","Percent_mCG","Num_CG_reads","mCG_reads","Weighted_mCG","Num_CH","mCH","Percent_mCH","Num_CH_reads","mCH_reads","Weighted_mCH"]) + '\n')
        elif c[1] != name:
            out.write('\t'.join([str(name),str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(CGr),str(mCGr),str(numpy.float64(mCGr)/numpy.float64(CGr)),str(CH),str(mCH),str(numpy.float64(mCH)/numpy.float64(CH)),str(CHr),str(mCHr),str(numpy.float64(mCHr)/numpy.float64(CHr))]) + '\n')
            name = c[1]
            CG = mCG = CGr = mCGr = CH = mCH = CHr = mCHHr = 0
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[5])
                    CHr = CHr + int(c[4])
                    mCHr = mCHr + int(c[3])
        elif c[1] == name:
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[5])
                    CHr = CHr + int(c[4])
                    mCHr = mCHr + int(c[3])
    out.write('\t'.join([str(name),str(CG),str(mCG),str(numpy.float64(mCG)/numpy.float64(CG)),str(CGr),str(mCGr),str(numpy.float64(mCGr)/numpy.float64(CGr)),str(CH),str(mCH),str(numpy.float64(mCH)/numpy.float64(CH)),str(CHr),str(mCHr),str(numpy.float64(mCHr)/numpy.float64(CHr))]) + '\n')
    out.close()

#Calculates average weighted methylation in plant specific contexts across bins for all input features for use in making
#metaplots
#Use: feature_metaplot_plant( "input allc file", "input bed file of features", "input required genome file for bedtools",
#"output file", ignore standedness (default = False), number of windows (default 60: 20 upstream, 20 downstream, 20 across
#feature), basepairs upstream and downstream (default is 2000 bp), cutoff on min number reads per site (default is 0)
#outputs tab-delimited file: column 1 = Bin number, column 2 = CG methylation level, column 3 = CHG methylation level, 
#column 4 = CHH methylation level
def feature_metaplot_plant(allc,features_bed,genome_file,metaplot_out,ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0):
    windows2 = int(windows/3)
    counter = 1
    bed_list = ["u","f","d"]
    mC_bed = allc2bed(allc)
    a = pandas.read_table(features_bed,header=None)
    if ignoreStrand is True:
        p_bed = pybedtools.BedTool.from_dataframe(a)
    else:
        p_bed = pybedtools.BedTool.from_dataframe(a[a[5]=='+'])
        n_bed = pybedtools.BedTool.from_dataframe(a[a[5]=='-'])
    out = open(metaplot_out, "w")
    CG = mCG = CHG = mCHG = CHH = mCHH = 0
    out.write('\t'.join(["Bin","mCG","mCHG","mCHH"]) + '\n')
    for y in bed_list:
        if y == "u":
            if ignoreStrand is True:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pybedtools.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pybedtools.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pybedtools.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        elif y == "f":
            if ignoreStrand is True:
                w_bed = pybedtools.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
            else:
                pw_bed = pybedtools.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
        elif y == "d":
            if ignoreStrand is True:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pybedtools.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pybedtools.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pybedtools.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        if ignoreStrand is True:
            w_bed = pybedtools.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pybedtools.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pandas.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,windows2+1)):
            for c in m.itertuples():
                if c[4].endswith("_"+str(x)):
                    if int(c[3]) >= int(cutoff):
                        if c[1].startswith("CN") or c[1].endswith("N"):
                            continue
                        elif c[1].startswith("CG"):
                            CG = CG + int(c[3])
                            mCG = mCG + int(c[2])
                        elif c[1].endswith("G"):
                            CHG = CHG + int(c[3])
                            mCHG = mCHG + int(c[2])
                        else:
                            CHH = CHH + int(c[3])
                            mCHH = mCHH + int(c[2])
            out.write('\t'.join([str(counter),str(numpy.float64(mCG)/numpy.float64(CG)),str(numpy.float64(mCHG)/numpy.float64(CHG)),str(numpy.float64(mCHH)/numpy.float64(CHH))]) + '\n')                
            counter = counter + 1
            CG = mCG = CHG = mCHG = CHH = mCHH = 0
    out.close()
    
#Calculates average weighted methylation in nonplant specific contexts across bins for all input features for use in making
#metaplots
#Use: feature_metaplot_nonplant( "input allc file", "input bed file of features", "input required genome file for bedtools",
#"output file", ignore standedness (default = False), number of windows (default 60: 20 upstream, 20 downstream, 20 across
#feature), basepairs upstream and downstream (default is 2000 bp), cutoff on min number reads per site (default is 0)
#outputs tab-delimited file: column 1 = Bin number, column 2 = CG methylation level, column 3 = CH methylation level
def feature_metaplot_nonplant(allc,features_bed,genome_file,metaplot_out,ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0):
    windows2 = int(windows/3)
    counter = 1
    bed_list = ["u","f","d"]
    mC_bed = allc2bed(allc)
    a = pandas.read_table(features_bed,header=None)
    if ignoreStrand is True:
        p_bed = pybedtools.BedTool.from_dataframe(a)
    else:
        p_bed = pybedtools.BedTool.from_dataframe(a[a[5]=='+'])
        n_bed = pybedtools.BedTool.from_dataframe(a[a[5]=='-'])
    out = open(metaplot_out, "w")
    CG = mCG = CH = mCH = 0
    out.write('\t'.join(["Bin","mCG","mCHH"]) + '\n')
    for y in bed_list:
        if y == "u":
            if ignoreStrand is True:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pybedtools.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pybedtools.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pybedtools.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        elif y == "f":
            if ignoreStrand is True:
                w_bed = pybedtools.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
            else:
                pw_bed = pybedtools.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
        elif y == "d":
            if ignoreStrand is True:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pybedtools.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pybedtools.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pybedtools.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pybedtools.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        if ignoreStrand is True:
            w_bed = pybedtools.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pybedtools.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pandas.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,windows2+1)):
            for c in m.itertuples():
                if c[4].endswith("_"+str(x)):
                    if int(c[3]) >= int(cutoff):
                        if c[1].startswith("CN") or c[1].endswith("N"):
                            continue
                        elif c[1].startswith("CG"):
                            CG = CG + int(c[3])
                            mCG = mCG + int(c[2])
                        else:
                            CHH = CHH + int(c[3])
                            mCHH = mCHH + int(c[2])
            out.write('\t'.join([str(counter),str(numpy.float64(mCG)/numpy.float64(CG)),str(numpy.float64(mCH)/numpy.float64(CH))]) + '\n')                
            counter = counter + 1
            CG = mCG = CH = mCH = 0
    out.close()
    
        
if __name__ == "__main__":
    pass