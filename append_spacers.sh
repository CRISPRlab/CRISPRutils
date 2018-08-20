#!/bin/bash
#USAGE: 
	# append_spacers.sh 
	# appends all files in current directory with filename *_spacers.fa* into spacers_concat.fasta

SPACERS=()
while IFS=  read -r -d $'\0'; do
    SPACERS+=("$REPLY")
done < <(find . -maxdepth 1 -name '*_spacers.fa*' -print0)
echo "Total spacer files: ${#SPACERS[@]}"

if [ ! -f ./spacers_contcat.fasta ]; then
	rm spacers_concat.fasta
fi


for i in "${SPACERS[@]}"
do 
	
	echo ${i}
	
	# BASE_PREFIX=$(echo ${i} | cut -d. -f2)
	# echo ${BASE_PREFIX}

	#insert filename into spacer comments
	# perl -pi -e "s/>/>${BASE_PREFIX:2}_/g" ${i}

	cat ${i} >> spacers_concat.fasta

done

echo "Appending completed successfully. Data dumped into spacers_concat.fasta"
exit 0