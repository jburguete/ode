t1:5.0434363883750303473e-01;
t2:9.8306664283175104242e-01;
b20:4.9153330200210544211e-01;
b21:4.9153334082964560031e-01;
a20:2.5400016673006989580e-02;
a21:9.7459998332699301042e-01;
c20:0.0000000000000000000e+00;
c21:5.0434367867696622208e-01;
t3:1.0000000000000000000e+00;
b30:3.3046247891000984505e-01;
b31:3.3046250501413161748e-01;
b32:3.3907501607585853749e-01;
a30:3.2769056020421107392e-01;
a31:0.0000000000000000000e+00;
a32:6.7230943979578892608e-01;
c30:0.0000000000000000000e+00;
c31:5.0434367859596119295e-01;
c32:5.0434367867696622219e-01;
b30+b31+b32-1;
b31*t1+b32*t2-1/2;
b31*t1^2+b32*t2^2-1/3;
a20+a21-1;
a20*c20+a21*t1-b20;
a21*c21-b21;
a30+a31+a32-1;
a30*c30+a31*t1+a32*b20-b30;
a31*c31+a32*b21-b31;
a32*c32-b32;