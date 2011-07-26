for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20; do
  #statistics --com com/liver-seg0${i}.com --label training/liver-seg0${i}.nhdr
  grep location com/liver-seg0$i.com | awk '((NR % 2) == 0) { print $2, $3, $4 }' > com/liver-seg0$i.loc
done
