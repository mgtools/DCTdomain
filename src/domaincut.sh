#!/bin/bash
start=$(date +%s)
temp=$( realpath "$0" )
d="$(dirname "$0")"
while IFS= read -r line
do
  array=($(echo $line | tr " " "\n"))
  cef="CE/${array[0]}.ce"
  if test -f "$cef"; then
     result=`${d}/RecCut --input CE/${array[0]}.ce --name ${array[0]}`
     echo "$result"
  fi
done

end=$(date +%s)
echo "#time used: $(($end-$start)) seconds"
