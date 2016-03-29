import math

from random import randint

def shanks_discrete_log(h, g, p, q = None):
  ''' Find x such that h = (g ** x) % p using the baby-step giant-step algorithm 
  developed by Daniel Shanks
  
  q gives the range of values to check.  By default, it is q
  '''
  # print("Solving {} = {} ** x mod {} via Shanks".format(h,g,p))
  if h > p:
    h %= p
  if g > p:
    g %= p
  if q == None:
    q = p

  n = 1 + math.floor(math.sqrt(q))
  g_list = []
  
  # generate all the powers of g mod p
  for i in range(0,n+1):
    g_list.append( pow(g,i,p) )

  g_inv = modular_inverse(g, p)
  g_neg_n = pow(g_inv,n,p)
  gamma = h % p

  # generate and check each power of h^-i*n sequentially
  for i in range(0,n+1):
    for j in range(0,n+1):
      if g_list[j] == gamma:
        return i*n + j
    gamma = (gamma * g_neg_n) % p

  print("ERROR")
  return -1

def ph_discrete_log(h,g,p):
  ''' Find x such that h = (g ** x) % p using the Pohlig-Hellman algorithm
  '''
  # print("Solving {} = {} ** x mod {} via Pohlig-Hellman".format(h,g,p))
  # first, get the factors of p-1
  factors = factorize(p-1)

  # then, build lists of smaller g_i's and h_i's
  gs = []
  hs = []
  for i in range(0,len(factors)):
    power = (p-1)//factors[i]
    gs.append(pow(g,power,p))
    hs.append(pow(h,power,p))


  # then, solve each smaller DLP
  xs = []
  for i in range(0,len(factors)):
    if gs[i] == 1:
      if hs[i] == 1:
        xs.append(0)
      else:
        print("ERROR: NO SOLUTION EXISTS")
        return -1
    else:
      xs.append(shanks_discrete_log(hs[i],gs[i],p,factors[i]))

  # print(xs)

  congruences = list(zip(xs,factors))

  # finally, use CRT to solve for x
  return chinese_remainder_theorem(congruences)
  

def modular_inverse(a, n):
  """ Find the modular inverse of a mod n 
  Based on pseudocode at 
  https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Modular_integers
  """
  t = 0
  new_t = 1
  r = n
  new_r = a

  while new_r != 0:
    quotient = r // new_r

    temp_t = t
    t = new_t
    new_t = temp_t - quotient * new_t

    temp_r = r
    r = new_r
    new_r = temp_r - quotient * new_r

  if r > 1:
    return -1
  elif t < 0:
    t = t + n

  return t

def gcd(a,b):
  ''' find the greatest common factor of a and b with the Euclidean Algorithm
  '''
  while b != 0:
    temp = b
    b = a % b
    a = temp

  return a

def isPrime(n,k=10):
  ''' determine if n is prime using the Miller-Rabin test
  k is the accuracy desired
  '''
  if n in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]:
    return True

  # write n - 1 = 2^r * d, where d is odd
  d = n-1
  r = 0
  while d % 2 == 0:
    r += 1
    d //= 2
  
  for i in range(0,k):
    a = randint(2,n-2)
    x = pow(a,d,n) # a**d % n

    if x==1 or x==n-1:
      continue

    for j in range(0,r-1):
      print(x)
      x = (x*x) % n

      if x==1:
        return False
      elif x==n-1:
        break

    if x==n-1:
      continue

    return False

  return True

def pollard_g(x,n):
  ''' helper function for pollard rho method
  '''
  return (x * x + 1) % n

def pollard_rho(n,x=2):
  ''' Finds the next factor of n using the Pollard rho method
  '''
  y=x
  d = 1

  while d==1:
    x = pollard_g(x,n)
    y = pollard_g(pollard_g(y,n),n)
    d = gcd(abs(x-y),n)

  if d == n:
    return -1
  else:
    return d

  
def factorize(n,x=2):
  ''' Return a list of tuples (factor, power) such that
  n = sum( factor**power) over the whole list
  '''
  factors = []
  while n > 1:
    q = pollard_rho(n,x)
    if q >= 2 and not isPrime(q):
      # print("refactoring {}".format(q))
      # when you refactor, try a different x value
      sub_factors = factorize(q,randint(0,q-1))
      factors = factors + sub_factors
    elif q == -1:
      q = n
      factors.append(q)
    else:
      factors.append(q)
    n = n // q

  return factors

def chinese_remainder_theorem(congruences):
  ''' Perform the CRT on the list of 2-tuples representing congruences;
  each tuple (a,p) represents x = a mod p.
  Returns a result unique modulo the product of all the p values
  '''
  big_p = 1
  for i in range(0,len(congruences)):
    big_p *= congruences[i][1]

  x = 0
  for i in range(0,len(congruences)):
    N = big_p // congruences[i][1]
    x += congruences[i][0] * N * modular_inverse(N,congruences[i][1])

  return x


