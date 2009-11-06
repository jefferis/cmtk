for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20; do
  statistics --com liver-seg0${i}-30.com --label ../resamp30/liver-seg0${i}-30.nhdr
  grep location liver-seg0$i-30.com | awk '((NR % 2) == 0) { print $3, $4, $5 }' > liver-seg0$i-30.loc
done
