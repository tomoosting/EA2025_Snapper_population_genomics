#!/bin/bash

lib_file=/nfs/scratch/oostinto/projects/snapper/resources/paleomix_lists/snapper_list.txt
paleomix_dir=/nfs/scratch/oostinto/projects/snapper/data/paleomix
blank_yaml=snapper_blank_1.0.yaml

N=$( wc -l < $lib_file )

while read -r line
do
    #get info
	sample=$( echo $line | cut -d ' ' -f 1 )
	extrac=$( echo $line | cut -d ' ' -f 2 )
	
	#make directories
	mkdir $paleomix_dir/$sample
	
	#copy blank yaml
	yaml=$paleomix_dir/$sample/$sample'_1.0.yaml'
	cp $blank_yaml $yaml
	
	#insert sample info
	sed -i "s/xxsamplexx/$sample/g" $yaml
	sed -i "s/xxextxx/$extrac/g" $yaml
	
done < $lib_file

