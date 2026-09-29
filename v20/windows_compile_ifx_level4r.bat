@REM  + + +
@echo off

setlocal

set COMPILE_FLAGS=/O2 /check:format /Qipo /Qprec-div- /QaxSSE4.2 /trace  /Qdiag-disable:8291 /extend-source

set LINK_FLAGS=/O2 /Qipo /check:format /Qprec-div- /QaxSSE4.2 


ifx /compile-only %COMPILE_FLAGS% level4r.f

ifx /exe:level4r.exe %LINK_FLAGS% level4r.obj


del *.obj
del *.mod
