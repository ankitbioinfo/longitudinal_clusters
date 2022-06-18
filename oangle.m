    
 function  angle= oangle(u,v)
      factor=180/pi;
      angle=factor*atan2(norm(cross(u,v)),dot(u,v));
 end