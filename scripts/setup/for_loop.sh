#!/usr/bin/env bash

WORDLIST="here is a list of words resulting in a sentence"

for word in $WORDLIST
do
  echo "'$word' with e's removed looks like:"
  echo $word | tr -d "e"
done
