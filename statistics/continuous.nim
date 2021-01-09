import ../math/special,
       math

## TODO: Add F-distribution

type
  ## ##
  ## Normal Distribution
  ## Parameters: mu,sigma (float)
  ## ## 
  Normal* = object
    mu,sigma: float 

proc StandardNormal*(): Normal =
  ## ##
  ## Standard Normal Distribution
  ## Parameters: mu (0.0), sigma(1.0) (float)
  ## ##
  Normal(mu: 0.0,sigma: 1.0)

proc mean*(n: Normal): float = n.mu

proc median*(n: Normal): float = n.mu

proc variance*(n: Normal): float =
  math.pow(n.sigma,2.0)

proc std*(n: Normal): float = n.sigma

proc mode*(n: Normal): float = n.mu

proc skewness*(n: Normal): float = 0.0

proc kurtosis*(n: Normal): float = 0.0

proc pdf*(n: Normal, x: float): float =
  var z: float = (x - n.mu) / (2.0*n.sigma)
  (1.0 / (sqrt(2.0*PI)*n.sigma))*exp(-z*z)

proc cdf*(n: Normal, x: float): float =
  (1.0 / 2.0)*(1.0 + erf((x - n.mu)/(n.sigma*math.sqrt(2.0))))

type 
  ## ##
  ## Uniform Distribution
  ## Parameters: a,b (float)
  ## ## 
  Uniform* = object  
    a,b: float 

proc StandardUniform*(): Uniform =
  ## ##
  ## Standard Uniform Distribution
  ## Parameters: a (0.0), b (1.0) (float) 
  ## ##
  Uniform(a: 0.0,b: 1.0)

proc mean*(u: Uniform): float =
  (u.a + u.b) / 2.0

proc variance*(u: Uniform): float =
  pow((u.b - u.a),2.0) / 12.0

proc std*(u: Uniform): float =
  sqrt(variance(u))

proc median*(u: Uniform): float =
  mean(u)

proc skewness*(u: Uniform): float = 0.0

proc kurtosis*(u: Uniform): float = -6/5

proc pdf*(u: Uniform,x: float): float =
  if u.a <= x and x <= u.b:
    result = 1.0/(u.b - u.a)
  else:
    result = 0.0

proc cdf*(u: Uniform, x: float): float =
  if u.a < 0.0:
    result = 0.0
  elif u.a <= x and x <= u.b:
    result = (x - u.a) / (u.b - u.a)
  else:
    result = 1.0

type  
  ## ##
  ## Exponential Distribution
  ## Parameters: lambda (float) 
  ## ##
  Exponential* = object  
    lambda: float 

proc StandardExponential*(): Exponential =
  ## ##
  ## Standard Exponential Distribution
  ## Parameters: lambda (1.0) (float) 
  ## ##
  Exponential(lambda: 1.0)

proc mean*(e: Exponential): float =
  1.0/e.lambda

proc variance*(e: Exponential): float =
  1.0/(e.lambda*e.lambda)

proc std*(e: Exponential): float =
  sqrt(variance(e))

proc median*(e: Exponential): float =
  ln(2.0)/e.lambda

proc skewness*(e: Exponential): float = 2.0

proc kurtosis*(e: Exponential): float = 6.0

proc pdf*(e: Exponential,x: float): float =
  ## Scaled version
  (1.0 / e.lambda) * exp(-e.lambda*x)

proc cdf*(e: Exponential, x: float): float =
  ## Scaled version
  1.0 - exp(-x/e.lambda)

type  
  ## ##
  ## Gamma Distribution (scaled version)
  ## Parameters: alpha,beta (float)  
  ## ##
  Gamma* = object   
    alpha,beta: float  

proc mean*(g: Gamma): float =
  g.alpha * g.beta

proc variance*(g: Gamma): float =
  g.alpha * g.beta * g.beta

proc mode*(g: Gamma): float =
  (g.alpha - 1.0) * g.beta

proc skewness*(g: Gamma): float =
  2.0 / sqrt(g.alpha)

proc kurtosis*(g: Gamma): float =
  6.0 / g.alpha

proc pdf*(g: Gamma,x: float): float =
  (1.0/(gamma(g.alpha)*pow(g.beta,g.alpha)))*pow(x,g.alpha - 1.0) * exp(-x / g.beta)

proc cdf*(g: Gamma, x: float): float =
  (1/(gamma(g.alpha)))*ligamma(g.alpha,x / g.beta)

type
  ## ##
  ## Beta Distribution
  ## Parameters: alpha,beta (float)
  ## ##
  Beta* = object
    alpha,beta: float

proc mean*(b: Beta): float =
  b.alpha / (b.alpha + b.beta)

proc variance*(b: Beta): float =
  var apb: float = b.alpha + b.beta ## alpha + beta term 
  (b.alpha * b.beta) / ((apb*apb)*(apb + 1.0))

proc pdf*(b: Beta,x: float): float =
  var 
    t0: float = (b.alpha-1)*ln(x) + (b.beta-1.0)*ln(1.0-x) - lbeta(b.alpha,b.beta)
  result = exp(t0)

proc cdf*(b: Beta,x: float): float =
  ## ##
  ## Regularized version
  ## ##
  incbeta(b.alpha,b.beta,x) / beta(b.alpha,b.beta)

type 
  ## ##
  ## Student T Distribution
  ## Parameters: nu (float)
  ## ##
  Student* = object
    nu: float

proc mean*(s: Student): float =
  if s.nu > 1.0:
    result = NaN
  result = 0.0

proc median*(s: Student): float = 0.0

proc mode*(s: Student): float = 0.0

proc variance*(s: Student): float =
  if s.nu > 2.0:
    result = s.nu / (s.nu - 2.0)
  result = Inf

proc skewness*(s: Student): float =
  if s.nu > 3.0:
    result = 0.0
  result = NaN
  
proc kurtosis*(s: Student): float =
  if 2.0 < s.nu and s.nu <= 4.0:
    result = Inf
  
  elif s.nu > 4.0:
    result = 6.0/(s.nu-4.0)
  
  else:
    result = NaN
  
proc pdf*(s: Student,x: float): float =
  var
    p: float = 0.5*(s.nu + 1.0)
    num: float = lgamma(p)
    den: float = 0.5*ln(s.nu*PI) + lgamma(0.5*s.nu)
    t0: float = -p*ln(1.0 + x*x/s.nu)
    
  result = exp(num - den + t0)
  
proc cdf*(s: Student,x: float): float =
  var
    p: float = 0.5*(s.nu + 1.0)
    q: float = 0.5*s.nu
    tmp_x: float = x
  if tmp_x < 0:
    tmp_x = -tmp_x
  var 
    t0: float = ln(tmp_x) + lgamma(p)
    t1: float = ln(hyp2f1(0.5,p,1.5,-tmp_x*tmp_x/s.nu))
    t2: float = 0.5*ln(PI*s.nu) + lgamma(q)
    v: float = t0 + t1 - t2
   
  result = 0.5 + exp(v)

type
  ## ##
  ## Weibull Distribution (2-parameter)
  ## Parameters: kappa (float) lambda (float)
  ## About: https://en.wikipedia.org/wiki/Weibull_distribution
  ## ##
  Weibull* = object
    kappa,lambda: float 

proc mean*(w: Weibull): float =
  w.lambda*gamma(1.0 + 1.0/w.kappa)

proc median*(w: Weibull): float =
  w.lambda*pow(ln(2.0),1.0/w.kappa)

proc mode*(w: Weibull): float =
  if w.lambda > 1.0:
    result = w.lambda*pow(((w.kappa - 1.0)/w.kappa),1/w.kappa)
  
  result = 0.0

proc variance*(w: Weibull): float =
  var 
    t0: float = gamma(1.0+2.0/w.kappa)
    t1: float = gamma(1.0+1.0/w.kappa)
  
  result = w.lambda*w.lambda*(t0 - t1*t1)
  
proc skewness*(w: Weibull): float =
  var 
    mu: float = w.mean()
    v: float = w.variance()
    t0: float = gamma(1.0+3.0/w.kappa)
    t1: float = w.lambda*w.lambda*w.lambda
    t2: float = 3.0*mu*v
    t3: float = mu*mu*mu
    t4: float = v*v*v
  
  result = (t0*t1 - t2 - t3)/t4

proc kurtosis*(w: Weibull): float =
  var 
    g1: float = gamma(1.0+1.0/w.kappa)
    g2: float = gamma(1.0+2.0/w.kappa)
    g3: float = gamma(1.0+3.0/w.kappa)
    g4: float = gamma(1.0+4.0/w.kappa)
    
    t0: float = g1*g1*g1*g1
    t1: float = g1*g1*g2
    t2: float = g2*g2
    t3: float = g1*g3
    t4: float = g4
    t5: float = g2-g1*g1
  
  result = (-6.0*t0+12.0*t1-3.0*t2-4.0*t3+t4) / t5*t5
    

