#tsc --outFile ff.js Atom.ts Bond.ts Angle.ts DihedralRB.ts DihedralCosine.ts \
#                    ForceField.ts System.ts VectorTools.ts main.ts

rm -f dat/tmp.xyz 
#tsc --outFile ff.js *.ts

#tsc --target ES3 --outFile ff.js --project .
tsc -p .
