


def E(n): return (3*n-8-(-2)**(3-n))/18

def R(n):
  m = (n-3)//2
  a = (1/2)*(1-1/2**m)
  b = 1/36*(
     +(3*n-11)*(2-1/2**m)
     -6*((-m+2**(m+1)-2)/2**m)
     +sum(2**(3-n+j+1) for j in range(m+1))
  )
  b = 1/36*(
     +(3*n-11)*(2-1/2**m)
     -6*((-m+2**(m+1)-2)/2**m)
     +(2**(m+1)-1)*2**(4-n)
  )
  return a+b