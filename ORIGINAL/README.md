To preprocess these files and remove my experimental code
(to make the code work like in the actual package), install
'rpl' and 'filepp'; then:

```
rpl '#ifndef' '#\!ifndef' *.R
rpl '#ifdef' '#\!ifdef' *.R
rpl '#else' '#\!else' *.R
rpl '#endif' '#\!endif' *.R
filepp -ov -kc '#!' -DPACKAGE *.R
```

