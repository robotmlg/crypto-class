#!/usr/bin/python3

import math

# import my CryptoMath library
from CryptoMath import *

def main():
  prime = 7975853687392939527669852234831196561598260079
  p_root = 23
  k_pub = 3063041257316300739044850830867868510529081573
  c_1 = 4119045320833673114078985432258112630764384218
  c_2 = 6257275893039194294649737427079765917239570228

  
  private_keys = elgamal_crack(prime, p_root, k_pub, c_1, c_2)
  
  print("Alice's private key: ",private_keys['a'])
  print("Bob's random key: ",private_keys['k'])
  print("Message: ",private_keys['m'])

def elgamal_crack(p, g, k_pub, c_1, c_2):
  """ Decode the elgamal cipher based on the given information, that is find
  a and k used to encode c_1 and c_2 (see below)

  p -- the prime base
  g -- a primative root of p
  k_pub -- equal to (g ** a) % p
  c_1 -- (g ** k) % p
  c_2 -- m * (k_pub ** k) % p

  where a is the private key of Alice, and  2 < k <= p-2 is a random number
  chosen by Bob
  """
  ret = {}  # a dictionary to hold the return values
  ret['a'] = 0
  ret['k'] = 0
  ret['m'] = 0

  # find the private keys with discrete logs
  k = ph_discrete_log(c_1,g,p)
  a = ph_discrete_log(k_pub,g,p)

  # decrypt the message
  x = pow(modular_inverse(c_1,p),a,p)
  m = (x * c_2) % p

  ret['a'] = a
  ret['k'] = k
  ret['m'] = decode_text(m)
  return ret

def decode_text(m):
  ''' convert the message to ASCII with the encoding scheme:
  a -> 10, b -> 11, c -> 12, etc, then concat the numbers, so
  abc -> 101112
  '''
  m_str = ""
  while m > 0:
    next_char_val = m % 100
    next_char = chr(next_char_val - 10 + ord('a'))
    m_str += next_char
    m //= 100

  return m_str[::-1]
    

if __name__ == "__main__":
  main()
