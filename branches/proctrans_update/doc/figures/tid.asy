import settings;
outformat="pdf";
size(15cm,0);
int[][] tid= {
  {1,1,1,0,0,2,2,2,0,0,0,3,3,0,0,0},
  {0,0,0,1,1,0,0,0,2,2,2,0,0,0,0,0},
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}
};
pen[][] tidfills= {
  {white, palered, mediumred, heavyred},
  {white, palegreen, mediumgreen, heavygreen},
  {white, paleblue, mediumblue, heavyblue}\
};
int nrows=tid.length;
int ncols=tid[0].length;
path sq=(0,0)--(0,1)--(1,1)--(1,0)--cycle;
for (int cr=0; cr<nrows; ++cr){
  label(string(cr+1),(0,nrows-cr-.5),W);
  for (int cc=0; cc<ncols; ++cc){
    filldraw(shift((cc,nrows-cr-1))*sq,tidfills[cr][tid[cr][cc]]);
    label(string(tid[cr][cc]),(cc+.5,nrows-cr-.5));
  }
}
for (int cc=0; cc<ncols; ++cc){
  label(string(cc+1),(cc+0.5,nrows),N);
}
label("Ensemble number",(ncols/2,nrows+.5),N);
label(rotate(90)*Label("Cross Section"),(-0.5,nrows/2),W);

draw((0,-5){right}..(4,-5)..{up}(5,-2)^^(7,-2){down}..{right}(10,-5)^^(0,-7){right}..(6,-7.5)..{right}(10,-7));
draw((3,-5.6)--(3,-7),palered+linewidth(2));
draw((3.2,-5.4)--(3.2,-7.2),mediumred+linewidth(2));
draw((2.8,-5.5)--(2.9,-7),heavyred+linewidth(2));
label("1",(3,-5.3),N);


draw((9.45,-5.3)--(9.5,-6.8),paleblue+linewidth(2));
label("3",(9.45,-5),N);


draw((5.1,-3)--(7.1,-3.1),palegreen+linewidth(2));
draw((5,-3.3)--(7.1,-3.3),mediumgreen+linewidth(2));
label("2",(5,-3),W);
