fpprec:20;
fpprintprec:fpprec;
t1:9.9999999999999999973b-01;
t2:5.0000000000000000022b-01;
b20:2.5000000000000000008b-01;
b21:2.5000000000000000014b-01;
a20:7.4999999999999999978b-01;
a21:2.5000000000000000022b-01;
c20:0.0000000000000000000b+00;
c21:9.9999999999999999973b-01;
t3:1.0000000000000000000b+00;
b30:1.6666666666666666685b-01;
b31:1.6666666666666666667b-01;
b32:6.6666666666666666647b-01;
a30:3.3333333333333333332b-01;
a31:0.0000000000000000000b+00;
a32:6.6666666666666666668b-01;
c30:3.6591823321385774650b-19;
c31:9.9999999938772709542b-01;
c32:9.9999999999999999967b-01;
b30+b31+b32+-1;
b31*t1+b32*t2+-1/2;
b31*t1^2+b32*t2^2+-1/3;
b32*(b21*t1+0)+-1/6;
b31*t1^3+b32*t2^3+-1/4;
a20+a21+-1;
a20*c20+a21*t1+-b20;
a21*c21+-b21;
a30+a31+a32+-1;
a30*c30+a31*t1+a32*b20+-b30;
a31*c31+a32*b21+-b31;
a32*c32+-b32;
