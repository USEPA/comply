@REM  + + +
@echo off

setlocal

set COMPILE_FLAGS=/O2 /check:format /Qipo /Qprec-div- /QaxSSE4.2 /trace  /Qdiag-disable:8291 /extend-source

set LINK_FLAGS=/O2 /Qipo /check:format /Qprec-div- /QaxSSE4.2 


ifx /compile-only %COMPILE_FLAGS% level4.f
ifx /compile-only %COMPILE_FLAGS% skrabel75.f

ifx /exe:level4.exe %LINK_FLAGS% level4.obj skrabel75.obj


del *.obj
del *.mod
