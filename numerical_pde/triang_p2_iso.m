function [c4n,n4e,Db] = triang_p2_iso()
phi = (0:pi/8:3*pi/2)';
c4n = [0,0;[cos(phi),sin(phi)]];
n4e = [1 2 4 0 3 0;1 4 6 0 5 0;1 6 8 0 7 0;1 8 10 0 9 0;
    1 10 12 0 11 0;1 12 14 0 13 0];
Db = [1 2 0;2 4 3;4 6 5;6 8 7;8 10 9;10 12 11;12 14 13];
