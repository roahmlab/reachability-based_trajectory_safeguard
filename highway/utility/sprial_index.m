function [x, y]= sprial_index(n)
  k = ceil((sqrt(n) - 1)/2);
  t = 2*k;
  m = (t + 1)^2;
  if n >= (m - t)
    x = k-(m-n); y = -k;
    return;
  else
    m = m -t;
  end
  
  if (n >= m - t)
    x = -k; y = -k + (m-n);
    return 
  else
    m = m - t;
  end

  if (n >= m - t)
    x = -k + (m-n); y = k;
    return;
  else
    x = k;y =  k-(m-n-t);
    return;
  end

end