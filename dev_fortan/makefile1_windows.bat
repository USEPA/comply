@REM                                                                    + + +
REM @echo off

setlocal

REM set COMPILE_FLAGS=/O2 /check:format /Qipo /Qprec-div- /QaxSSE4.2 /trace  /Qdiag-disable:8291
set COMPILE_FLAGS=/extend-source 
REM set LINK_FLAGS=/O2 /Qipo /check:format /Qprec-div- /QaxSSE4.2 

ifx /compile-only %COMPILE_FLAGS% level1.f  

ifx /exe:level1.exe %LINK_FLAGS% level1.obj

del *.obj
del *.mod






