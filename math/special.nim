import math
import rationals

const
  EM = 0.5772156649015328606651209008240243 #Euler Constant
  P = 1.6449340668482264365 #Trigamma Constant 
  EPS = 1e-15

proc hyp1f1*(a,b,z: float): float =
  ## ## 
  ## Hypergeometric Function (1F1)
  ## Parameters: a,b,z (float)
  ## Source: Pearson (2008)
  ## Computation of Hypergeometric Functions
  ## ##
  const
    LIMIT = 100
  var 
    a0: float = 0.0
    b0,g0,u0: float = 1.0
    aj,bj,gj,uj: float
  
  for j in 1..LIMIT:
    aj = (a0 + b0)*j.float*(b+j.float-1.0)
    bj = b0*(a+j.float-1.0)*z
    gj = g0*j.float*(b+j.float-1.0)
    uj = (aj + bj) / gj
    
    if abs(u0 - uj) < EPS: break

    a0 = aj
    b0 = bj
    g0 = gj
    u0 = uj
  
  result = uj

proc hyp2f1*(a,b,c,z: float): float =
  ## ##
  ## Hypergeometric Function (2F1)
  ## Parameters: a,b,c,z (float)
  ## Source: Pearson (2008)
  ## Computation of Hypergeometric Functions
  ## ## 
  const 
    LIMIT = 100
  var 
    a0: float = 0.0
    b0,g0,u0: float = 1.0 
    aj,bj,gj,uj: float 

  for j in 1..LIMIT:
    aj = (a0 + b0) * j.float * (c+j.float-1.0)
    bj = b0 * (a+j.float-1.0)*(b+j.float-1.0)*z
    gj = g0 * j.float * (c + j.float - 1.0)
    uj = (aj + bj) / gj

    if abs(u0 - uj) < EPS: break
   
    a0 = aj
    b0 = bj
    g0 = gj
    u0 = uj
  
  result = uj

proc ligamma*(a, z: float): float =
  ## ##
  ## Lower Incomplete Gamma Function
  ## Source: https://wikipedia.org/wiki/Gamma_distribution
  ## ##
  var 
    t0: float = hyp1f1(a,1.0+a,-z)
    t1: float = a*ln(z)-ln(a)
    t2: float = t0+ln(t1)
  
  result = exp(t2)

proc uigamma*(a, z: float): float =
  ## ##
  ## Upper Incomplete Gamma Function
  ## Parameters: a,z (float)
  ## Source: https://wikipedia.org/wiki/Gamma_distribution
  ## This is based on the relation: 
  ## gamma(z) = ligamma(a,z) + uigamma(a,z)
  ## ##
  gamma(a) - ligamma(a,z)

proc digamma*(z: float): float =
  ## ##
  ## Digamma Function approximation 
  ## Source: Jose Bernardo Algorithm AS 103 (1976) 
  ## Psi(Digamma) Function, Applied Statistics, Volume 25,
  ## Number 3 
  ## ## 
  const
    s: float = 1e-6
    s3: float = 1.0 / 12.0
    s4: float = 1.0 / 120.0
    s5: float = 1.0 / 252.0
    s6: float = 1.0 / 240.0
    s7: float = 1.0 / 132.0
  var
    r: float = 0.0
    tmp_z: float = z
  
  if z < 0:
    r = digamma(1.0-tmp_z) + PI/tan(-PI*tmp_z)
  if z <= s:
    r = -EM - (1.0/tmp_z) + (P*tmp_z)
  while tmp_z < 12.0:
    r -= 1.0/tmp_z
    tmp_z+=1
  if tmp_z >= 12.0:
    var 
      k: float = 1.0/tmp_z
    r += ln(tmp_z) - (0.5*k)
    k *= k
    r -= k*(s3-(k*(s4-(k*(s5-(k*(s6-(k*s7))))))))
  
  result = r

proc lbeta*(x, y: float): float =
  ## ##
  ## Log Beta Function
  ## ##
  lgamma(x) + lgamma(y) - lgamma(x+y)

proc beta*(x, y: float): float =
  ## ##
  ## Beta Function
  ## ##
  exp(lbeta(x,y))

proc incbeta*(a,b,x: float): float =
  ## ##
  ## Incomplete Beta Function
  ## Parameters: a,b (float)
  ## Taking advantage of the incomplete beta hypergeometric representation
  ## Source: https://dlmf.nist.gov/8.17
  ## ##
  var 
    t0: float = a*ln(x)
    t1: float = ln(a)
    t2: float = hyp2f1(a,1.0-b,a+1.0,x)
    v: float = t0 - t1 + ln(t2)
  
  result = exp(v)


