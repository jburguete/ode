fpprec:20;
fpprintprec:fpprec;
t1:4.9999999038322351139b-01;
t2:5.0000000961677685857b-01;
b20:-1.9233559465524265736b-08;
b21:5.0000002885033632409b-01;
t3:1.0000000000000000000b+00;
b30:3.8467106324379350679b-08;
b31:1.0016944064254999796b-14;
b32:9.9999996153288365868b-01;
t4:1.0000000000000000000b+00;
b40:1.6666666666666654341b-01;
b41:3.3333333974287389376b-01;
b42:3.3333332692379301948b-01;
b43:1.6666666666666654333b-01;
b40+b41+b42+b43+-1;
b41*t1+b42*t2+b43*t3+-1/2;
b41*t1^2+b42*t2^2+b43*t3^2+-1/3;
b42*(b21*t1+0)+b43*(b31*t1+b32*t2+0)+-1/6;
b41*t1^3+b42*t2^3+b43*t3^3+-1/4;
b43*(b32*(b21*t1+0)+0)+-1/24;
b42*(b21*t1^2+0)+b43*(b31*t1^2+b32*t2^2+0)+-1/12;
b42*t2*(b21*t1+0)+b43*t3*(b31*t1+b32*t2+0)+-1/8;
b41*t1^4+b42*t2^4+b43*t3^4+-1/5;
