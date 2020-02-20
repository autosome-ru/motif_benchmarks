#!/usr/bin/env sh
/app/fasta_multiline_to_singleline.sh $1 \
  | sed -re 's/\t/%tab%/g' \
  | awk 'NR%2{printf("%s\t",$0);next;} 1{print$0}' \
  | cat -n | sed -re 's/^\s+//' \
  | awk -F $'\t' '!_[$3]++' \
  | cut -f2- -d $'\t' \
  | awk -F $'\t' -e '{ print($1 "\n" $2) }' \
  | sed -re 's/%tab%/\t/g'
