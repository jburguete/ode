t1:4.9034468279033899312e-01;
t2:4.9116787655614316382e-01;
b20:7.9429357126352959691e-05;
b21:4.9108844719901681086e-01;
a20:6.0289829261058734922e-01;
a21:3.9710170738941265078e-01;
c20:-3.2283601359268342606e-01;
c21:1.2366817821748554707e+00;
t3:6.9220111079365015031e-01;
b30:1.6299315987005829812e-01;
b31:1.4671285852566361427e-01;
b32:3.8249509239792823791e-01;
a30:7.3891601369247352751e-01;
a31:6.6550961833103126278e-02;
a32:1.9453302447442334625e-01;
c30:1.7640001779929919764e-01;
c31:7.6903377799317240715e-01;
c32:1.9662218969315290451e+00;
t4:1.0000000000000000000e+00;
b40:1.4129167618693456281e-02;
b41:1.3250801947010741576e-02;
b42:9.0477673140511202213e-01;
b43:6.7843299029183780023e-02;
a40:6.2525768271313138388e-01;
a41:-1.0534960183092158471e-01;
a42:4.4570627641513473706e-01;
a43:3.4385642702655463827e-02;
c40:9.6195152258351580645e-02;
c41:1.9997723153792321682e+00;
c42:2.0004752883307635786e+00;
c43:1.9730123882181942463e+00;
t5:1.0000000000000000000e+00;
b50:1.6492271910621553779e-01;
b51:3.1966127307814364163e-01;
b52:3.1662582641040980437e-01;
b53:3.5902734419927496679e-02;
b54:1.6288744698530351954e-01;
a50:4.1819476343641609576e-01;
a51:1.1122251396547420700e-01;
a52:1.4708503388987956998e-01;
a53:2.2196434753801558282e-01;
a54:1.0153334117021454442e-01;
c50:1.7398691300937851156e-01;
c51:1.9197470229119062424e+00;
c52:9.5088225184008479118e-01;
c53:1.3071638718247413635e-01;
c54:1.6042754538357257843e+00;
a50+a51+a52+a53+a54-1;
b51*t1+b52*t2+b53*t3+b54*t4-1/2;
b52*b21*t1+b53*(b31*t1+b32*t2)+b54*(b41*t1+b42*t2+b43*t3)-1/6;
b51*t1^2+b52*t2^2+b53*t3^2+b54*t4^2-1/3;
b53*b32*b21*t1+b54*(b42*b21*t1+b43*(b31*t1+b32*t2))-1/24;
b52*b21*t1^2+b53*(b31*t1^2+b32*t2^2)+b54*(b41*t1^2+b42*t2^2+b43*t3^2)-1/12;
b52*t2*b21*t1+b53*t3*(b31*t1+b32*t2)+b54*t4*(b41*t1+b42*t2+b43*t3)-1/8;
b51*t1^3+b52*t2^3+b53*t3^3+b54*t4^3-1/4;
