%moving average of lad
%

function lad_mov = movlad(lad, ww)
x = lagg(lad, ww); 
lad_mov = mean(x,2); 