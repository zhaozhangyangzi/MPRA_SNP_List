#!/bin/bash 

ExtractSnpLD(){
  if [ $# -eq 4 ]
  then
    cut -f $3 $1 | grep -n -w -Ff $2 | cut -f 1 -d ":" | cat | while read line
    do
    head -"$line" $1 | tail -1 >> $4
    done
  else
  echo "dataFile leadSnpFile columnNumber outputFile"
  fi
}
