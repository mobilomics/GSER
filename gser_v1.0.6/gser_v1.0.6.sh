#!/bin/bash

#GSER v1.0.6

#change path to R if installed in another custom directory
pathToR=$(which R)

gser_sourcedir=${BASH_SOURCE[0]}

inputfile=$1
threads=$2
outputdir=$3

kmerRange=$(head -n 1 $inputfile|awk '{print $2}')
kmerRangeForLoop=${kmerRange//,/ }

#outputdir="all"
outputfilename="default_test"

Rscript_filename="generate_plots.R"

if [ ! -d "$outputdir" ]; then
	echo "$outputdir doesn't exist. Will be created."
	mkdir "$outputdir"
fi
outputdir_fullpath=$(readlink -f $outputdir)

inputfile=$(readlink -f $inputfile)
echo ">GSER started with file $inputfile"

cd $outputdir_fullpath
fileGrouping_file="fileGrouping.txt"
tail -n +2 $inputfile > $fileGrouping_file

groups=$(awk 'BEGIN{FS="\t"}{print $1}' $fileGrouping_file|sort -u|tr '\n' ' ')
echo $groups


echo "library(ggplot2)" > $Rscript_filename
echo "options(scipen=9)" >> $Rscript_filename
#echo "pdf(file=\""$outputfilename"_genomesizeestimations.pdf\",12,7)" >> $Rscript_filename

echo "overview_table <- data.frame()" >> $Rscript_filename
echo "plots <- list()" >> $Rscript_filename
echo "plot_count <- 1" >> $Rscript_filename

echo "find_peaks <- function(input_df,windowsize,max_kfrequency,rthreshold) { " >> $Rscript_filename
echo "        peak_counter <- 1 " >> $Rscript_filename
echo "        peaks <- numeric(0) " >> $Rscript_filename
echo "        for (x in 3:(max_kfrequency-windowsize)) { " >> $Rscript_filename
echo "          if(input_df\$V2[x+floor(windowsize/2)]>input_df\$V2[x] & input_df\$V2[x+floor(windowsize/2)]>input_df\$V2[x+windowsize]) { " >> $Rscript_filename
echo "                  peak_fitting <- lm(input_df\$V2[(x):(x+windowsize)] ~ poly((x):(x+windowsize),2,raw=TRUE)) " >> $Rscript_filename
echo "                  peak_fitting_summary <- summary(peak_fitting) " >> $Rscript_filename
echo "                  if(peak_fitting_summary\$adj.r.squared>=rthreshold){ " >> $Rscript_filename
echo "                          peak <- which(input_df\$V2==max(input_df\$V2[(x):(x+windowsize)])) " >> $Rscript_filename
echo "                          if(length(grep(peak,peaks))==0){ " >> $Rscript_filename
echo "                                  print(paste(x,x+windowsize,peak_fitting_summary\$adj.r.squared)) " >> $Rscript_filename
echo "                                  print(peak) " >> $Rscript_filename
echo "                                  peaks[peak_counter] <- peak " >> $Rscript_filename
echo "                                  peak_counter <- peak_counter+1  " >> $Rscript_filename
echo "                          } " >> $Rscript_filename
echo "                  } " >> $Rscript_filename
echo "          } " >> $Rscript_filename
echo "  " >> $Rscript_filename
echo "  } " >> $Rscript_filename
echo "        selected_peak <- max(peaks) " >> $Rscript_filename
echo "        return(selected_peak) " >> $Rscript_filename
echo "} " >> $Rscript_filename

echo "Obtaining lengths for eachlibrary using sample size N=1000"
for file in $(awk '{print $2}' fileGrouping.txt)
do
	echo $file
	$pathToR/Rscript $gser_sourcedir/getsamplemean_v1.0.R $file 1000 >> temp_lengthinfo.txt 2>/dev/null
done

paste fileGrouping.txt temp_lengthinfo.txt > fileGrouping_lengthinfo.txt


echo "fileGrouping_lengthinfo <- read.delim(\"fileGrouping_lengthinfo.txt\",header=FALSE,sep=\"\\t\")" >> $Rscript_filename

for group in $groups
do
	group_filename="GSER_"$group
	group_filename=$outputdir_fullpath"/"$group_filename

	#TO DO
	#store each group filepaths into a file
	echo -e "awk -v ID=\"$group\" 'BEGIN{FS=OFS=\"\\t\"}(\$1==ID){print \$0}' $fileGrouping_file > $group_filename"
	awk -v ID="$group" 'BEGIN{FS=OFS="\t"}($1==ID){print $0}' $fileGrouping_file > $group_filename

	group_filename_list=$group_filename"_fastqList"
	echo -e "awk 'BEGIN{FS=OFS=\"\t\"}{print \$2}' $group_filename > $group_filename_list"
	awk 'BEGIN{FS=OFS="\t"}{print $2}' $group_filename > $group_filename_list

	#run ntcard
	outputbasename="GSER_ntCard_"$group
	#outputbasename="GSER_ntcard_"$group
        ntcardlog=$outputbasename".ntcardlog"
        cmd=">CMD: /home/bvaldebenito/NTCARD/ntCard-master/ntcard -p $outputbasename -t$threads -k$kmerRange @$group_filename_list > $ntcardlog"
	echo $cmd
	ntcard -p $outputbasename -t$threads -k$kmerRange @$group_filename_list > $ntcardlog

	echo "groupLength <- sum(fileGrouping_lengthinfo[which(fileGrouping_lengthinfo\$V1==\"$group\"),]\$V4)" >> $Rscript_filename


	
	for k in $kmerRangeForLoop
        do
                histfilename=$outputbasename"_k"$k".hist"
                #basename=${outputbasename/GSER_ntCard_/}
		basename=$group"_k"$k
                basename=${basename//-/_}
                skip=2

		hist_f1_f0=$(head -n 2 $histfilename|awk 'BEGIN{sumf=0}{sumf=sumf+$2}END{print sumf}')
		if [[ $hist_f1_f0 == 0 ]]
		then
        		continue
		fi

		frequency_k1to100=$(tail -n +3 $histfilename|head -n 100|awk 'BEGIN{sumfreqs=0}{sumfreqs=sumfreqs+$2}END{print sumfreqs}')		
		if [[ $frequency_k1to100 == 0 ]]
                then
                        continue
                fi


		echo "$basename <- read.delim(\"$histfilename\",skip=$skip,header=FALSE,sep=\"\t\")"

		echo "input_df <- $basename"
		echo "#try1"
		echo "selected_peak <- find_peaks(input_df,25,1000,0.9)"
		echo "#try2"
		echo "if(selected_peak==-Inf){"
		echo "        selected_peak <- find_peaks(input_df,10,25,0.9)"
		echo "}"
		echo "#try3"
		echo "if(selected_peak==-Inf){"
                echo "        selected_peak <- find_peaks(input_df,10,25,0.85)"
                echo "}"
		echo "if(selected_peak==-Inf){"
                echo "        selected_peak <- find_peaks(input_df,10,25,0.8)"
                echo "}"
		echo "#no luck"
		echo "if(selected_peak==-Inf){"
		echo "        selected_peak <- find_peaks(input_df,5,15,0.8)"
		echo "}"
		echo "if(selected_peak<=100 & selected_peak!=-Inf){xmax <- 100}"
		echo "if(selected_peak>100 & selected_peak<=500){xmax <- 500}"

		echo "max_k <- selected_peak"
#                echo "$basename <- read.delim(\"$histfilename\",skip=$skip,header=FALSE,sep=\"\t\")"
#                echo "subset_to_query_valley <- $basename[1:10,]"
#                echo "min_kmer_ocurrence <- subset_to_query_valley\$V1[which(subset_to_query_valley\$V2==min(subset_to_query_valley\$V2))[1]]"
#                #echo "$basename <- $basename[2:1000,]"
#                echo "$basename <- $basename[1:100,]"
#                echo "max_k <- $basename\$V1[which($basename\$V2==max($basename\$V2))][1]"
#
		echo "if(selected_peak!=-Inf){"
		echo "estimated_kmer_cov <- selected_peak"	
                echo "estimated_kmer_number <- sum(as.double($basename\$V1)*as.double($basename\$V2))"
                echo "genome_size <- sum(as.double($basename\$V1)*as.double($basename\$V2))/max_k"
                echo "rounded_genome_size <- round(genome_size/1000000000,digits=5)"
                echo 'genome_size_label <- paste("Estimated genome size: ",rounded_genome_size," Gbp",sep="")'
                echo 'max_k_label <- paste("Estimated coverage: ",max_k,sep="")'
                echo 'estimated_kmer_number_label <- paste("Estimated K-mer number: ",estimated_kmer_number,sep="")'
                echo 'plot_label <- paste(genome_size_label,max_k_label,estimated_kmer_number_label,sep="\n")'
                #echo "plots[[plot_count]] <- ggplot($basename)+geom_line(aes(x=V1,y=V2,group=1))+theme_classic()+labs(y=\"Frequency\",x=\"K-mer ocurrence\",title=\"$basename\")+annotate(\"text\",label=plot_label,x=70,y=max($basename\$V2))+geom_vline(xintercept=estimated_kmer_cov,linetype=\"dashed\",color=\"red\")"
		echo "plots[[plot_count]] <- ggplot(input_df)+geom_line(aes(x=V1,y=V2,group=1))+theme_classic()+labs(y=\"Frequency\",x=\"K-mer ocurrence\",title=\"$basename\")+geom_vline(xintercept=selected_peak,linetype=\"dashed\",colour=\"red\")+scale_y_continuous(limits=c(0,input_df\$V2[selected_peak]*1.25))+scale_x_continuous(limits=c(0,xmax))+annotate(\"text\",label=plot_label,x=xmax*0.9,y=input_df\$V2[selected_peak]*1.2)"
		echo "current_stats <- data.frame(A=\"$group\",B=$k,C=estimated_kmer_cov,D=rounded_genome_size,E=groupLength)"
		echo "overview_table <- rbind(overview_table,current_stats)"
		echo "} else {"
                echo "	plots[[plot_count]] <- ggplot(input_df)+geom_line(aes(x=V1,y=V2,group=1))+theme_classic()+labs(y=\"Frequency\",x=\"K-mer ocurrence\",title=\"$basename\")"
		echo "}"

		echo "plot_count <- plot_count+1"
#		echo "print(\">$group $k finished ok\")"
                echo -e "\n\n"



        done >> $Rscript_filename



	
done


echo "library(reshape2)" >> $Rscript_filename
echo "colnames(overview_table) <- c(\"Group\",\"K\",\"Kcov\",\"GSE\")" >> $Rscript_filename
echo "reshaped <- melt(data=overview_table,id.vars=c(\"Group\",\"K\"),measure.vars=c(\"Kcov\",\"GSE\"))" >> $Rscript_filename
echo "reshaped\$variable <- factor(reshaped\$variable,levels=c(\"GSE\",\"Kcov\"))" >> $Rscript_filename


echo "width <- length(levels(as.factor(reshaped\$K)))*0.4*2" >> $Rscript_filename
echo "reshaped\$K <- as.factor(reshaped\$K)" >> $Rscript_filename
echo "pdf(file=\"gseR_overview.pdf\",width,7)" >> $Rscript_filename
echo "ggplot(reshaped,aes(x=K,y=value,colour=Group))+geom_point()+geom_line(aes(group=Group))+theme_classic()+labs(y=\"Estimated K-mer coverage\")+theme(legend.position=\"bottom\")+guides(col=guide_legend(nrow=3,byrow=TRUE))+facet_wrap(~variable,scales=\"free_y\",labeller=as_labeller(c(Kcov=\"Estimated K-mer Coverage\",GSE=\"Estimated Genome Size (Gbp)\")),strip.position=\"left\")+theme(strip.background=element_blank(),strip.placement=\"outside\",strip.text=element_text(size=11))+labs(x=\"K\",y=NULL)" >> $Rscript_filename
echo "dev.off()" >> $Rscript_filename


echo "colnames(overview_table) <- c(\"Group\",\"K\",\"Kcov\",\"GSE\",\"readLength\")" >> $Rscript_filename
echo "overview_table\$Gcov <- overview_table\$Kcov*overview_table\$readLength/(overview_table\$readLength-overview_table\$K+1)" >> $Rscript_filename
echo "reshaped2 <- melt(data=overview_table,id.vars=c(\"Group\",\"K\"),measure.vars=c(\"Gcov\",\"GSE\"))" >> $Rscript_filename
echo "reshaped2\$K <- as.factor(reshaped2\$K)" >> $Rscript_filename
echo "reshaped2\$variable <- factor(reshaped2\$variable,levels=c(\"GSE\",\"Gcov\"))" >> $Rscript_filename
echo "pdf(file=\"gseR_overview_gcov.pdf\",width,7)" >> $Rscript_filename
echo "ggplot(reshaped2,aes(x=K,y=value,colour=Group))+geom_point()+geom_line(aes(group=Group))+theme_classic()+labs(y=\"Estimated Genome coverage\")+theme(legend.position=\"bottom\")+guides(col=guide_legend(nrow=3,byrow=TRUE))+facet_wrap(~variable,scales=\"free_y\",labeller=as_labeller(c(Gcov=\"Estimated Genome coverage\",GSE=\"Estimated Genome Size (Gbp)\")),strip.position=\"left\")+theme(strip.background=element_blank(),strip.placement=\"outside\",strip.text=element_text(size=11))+labs(x=\"K\",y=NULL)" >> $Rscript_filename
echo "dev.off()" >> $Rscript_filename

#echo "ggplot(reshaped,aes(x=K,y=value,colour=Group))+geom_point()+geom_line()+theme_classic()+scale_colour_brewer(palette=\"Set1\")+labs(y=\"Estimated K-mer coverage\")+theme(legend.position=\"bottom\")+guides(col=guide_legend(nrow=3,byrow=TRUE))+facet_wrap(~variable,scales=\"free_y\",labeller=as_labeller(c(Kcov=\"Estimated K-mer Coverage\",GSE=\"Estimated Genome Size (Gbp)\")),strip.position=\"left\")+theme(strip.background=element_blank(),strip.placement=\"outside\",strip.text=element_text(size=11))+scale_x_continuous(breaks=c($kmerRange))+ylab(NULL)" >> $Rscript_filename


echo "pdf(file=\"gseR_kmer_histograms.pdf\",12,7)" >> $Rscript_filename

echo "for (i in 1:(plot_count-1)) {" >> $Rscript_filename
echo "	print(plots[[i]])" >> $Rscript_filename
echo "}" >> $Rscript_filename

echo "dev.off()" >> $Rscript_filename
echo "print(overview_table)" >> $Rscript_filename

echo "if(require(plotly)){" >> $Rscript_filename
echo "library(plotly)" >> $Rscript_filename
echo "ggplotly_obj <- ggplotly(ggplot(reshaped2,aes(x=K,y=value,colour=Group))+geom_point()+geom_line(aes(group=Group))+theme_classic()+labs(y=\"Estimated Genome coverage\")+theme(legend.position=\"bottom\")+guides(col=guide_legend(nrow=3,byrow=TRUE))+facet_wrap(~variable,scales=\"free_y\",labeller=as_labeller(c(Gcov=\"Estimated Genome coverage\",GSE=\"Estimated Genome Size (Gbp)\")),strip.position=\"left\")+theme(strip.background=element_blank(),strip.placement=\"outside\",strip.text=element_text(size=11))+labs(x=\"K\",y=NULL))" >> $Rscript_filename
echo "htmlwidgets::saveWidget(widget=ggplotly_obj,\"plotly_test.html\", selfcontained = FALSE)" >> $Rscript_filename
echo "}" >> $Rscript_filename


chmod a-rwx $Rscript_filename
chmod u+r $Rscript_filename
$pathToR/Rscript $Rscript_filename
rm -Rf $Rscript_filename
