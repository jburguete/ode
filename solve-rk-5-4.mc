solve([b51*t1+b52*t2+b53*t3+b54*t4=1/2],[b51]);
solve([b51*t1+b52*t2+b53*t3+b54*t4=1/2,
       b51*t1^2+b52*t2^2+b53*t3^2+b54*t4^2=1/3],
       [b52,b51]);
solve([b51*t1+b52*t2+b53*t3+b54*t4=1/2,
       b51*t1^2+b52*t2^2+b53*t3^2+b54*t4^2=1/3,
       b51*t1^3+b52*t2^3+b53*t3^3+b54*t4^3=1/4],
       [b53,b52,b51]);
solve([b51*t1+b52*t2+b53*t3+b54*t4=1/2,
       b51*t1^2+b52*t2^2+b53*t3^2+b54*t4^2=1/3,
       b51*t1^3+b52*t2^3+b53*t3^3+b54*t4^3=1/4,
       b53*b32*b21*t1+b54*(b42*b21*t1+b43*(b31*t1+b32*t2))=1/24],
       [b54,b53,b52,b51]);
solve([b51*t1+b52*t2+b53*t3+b54*t4=1/2,
       b51*t1^2+b52*t2^2+b53*t3^2+b54*t4^2=1/3,
       b51*t1^3+b52*t2^3+b53*t3^3+b54*t4^3=1/4,
       b53*b32*b21*t1+b54*(b42*b21*t1+b43*(b31*t1+b32*t2))=1/24,
       b52*b21*t1+b53*(b31*t1+b32*t2)+b54*(b41*t1+b42*t2+b43*t3)=1/6],
       [b41,b54,b53,b52,b51]);
solve([b51*t1+b52*t2+b53*t3+b54*t4=1/2,
       b52*b21*t1+b53*(b31*t1+b32*t2)+b54*(b41*t1+b42*t2+b43*t3)=1/6,
       b51*t1^2+b52*t2^2+b53*t3^2+b54*t4^2=1/3,
       b53*b32*b21*t1+b54*(b42*b21*t1+b43*(b31*t1+b32*t2))=1/24,
       b52*b21*t1^2+b53*(b31*t1^2+b32*t2^2)+b54*(b41*t1^2+b42*t2^2+b43*t3^2)
       =1/12,
       b52*t2*b21*t1+b53*t3*(b31*t1+b32*t2)+b54*t4*(b41*t1+b42*t2+b43*t3)=1/8,
       b51*t1^3+b52*t2^3+b53*t3^3+b54*t4^3=1/4],
       [b54,b53,b52,b51,b43,b42,b41]);
