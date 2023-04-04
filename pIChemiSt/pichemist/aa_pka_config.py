# module with hardcoded smarts of capped aminoacids 
D_cappedAA_smarts = {
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]OC(=O)C([H])([H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[O-]C(=O)C([H])([H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])O[H])N([H])C(=O)C([H])([H])[H]':   'S',
'[H]C([H])([H])C(=O)C([H])(C([H])(C([H])([H])[H])O[H])N([H])C(=O)C([H])([H])[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'N',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])S[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'A',
'[H]C([H])([H])C(=O)C([H])(C([H])(C([H])([H])[H])C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)C([H])([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)C([H])(C([S-])([H])[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])([H])N([H])C(=O)C([H])([H])[H]':   'G',
'[H]C([H])([H])C(=O)C1([H])C([H])([H])C([H])([H])C([H])([H])N1C(=O)C([H])([H])[H]':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'I',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]OC(=O)C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[O-]C(=O)C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'D',
'[H]N([H])C(=O)C([H])(C([H])([H])O[H])N([H])C(=O)C([H])([H])[H]':   'S',
'[H]C([H])([H])C([H])(O[H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'N',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]N([H])C(=O)C([H])(C([H])([H])S[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]N([H])C(=O)C([H])(C([H])([H])[H])N([H])C(=O)C([H])([H])[H]':   'A',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)N([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]N([H])C(=O)C([H])(C([S-])([H])[H])N([H])C(=O)C([H])([H])[H]':   'C',
'[H]N([H])C(=O)C([H])([H])N([H])C(=O)C([H])([H])[H]':   'G',
'[H]C([H])([H])C(=O)N1C([H])([H])C([H])([H])C([H])([H])C1([H])C(=O)N([H])[H]':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)N([H])[H])N([H])C(=O)C([H])([H])[H]':   'I'
}


D_NtermfreeAA_smarts = {
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])(C([H])([H])[H])O[H]':   'T',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C(=O)N([H])[H]':   'N',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])[H])N([H])[H]':   'A',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'V',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)C([H])(C([S-])([H])[H])N([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])([H])N([H])[H]':   'G',
'[H]C([H])([H])C(=O)C1([H])C([H])([H])C([H])([H])C([H])([H])N1[H]':   'P',
'[H]C([H])([H])C(=O)C([H])(N([H])[H])C([H])(C([H])([H])[H])C([H])([H])C([H])([H])[H]':   'I',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])(C([H])([H])[H])O[H]':   'T',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C(=O)N([H])[H]':   'N',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])C(=O)N([H])[H]':   'Q',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)C([H])(C([H])([H])[H])[N+]([H])([H])[H]':   'A',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'V',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)C([H])(C([S-])([H])[H])[N+]([H])([H])[H]':   'C',
'[H]C([H])([H])C(=O)C([H])([H])[N+]([H])([H])[H]':   'G',
'[H]C([H])([H])C(=O)C1([H])C([H])([H])C([H])([H])C([H])([H])[N+]1([H])[H]':   'P',
'[H]C([H])([H])C(=O)C([H])([N+]([H])([H])[H])C([H])(C([H])([H])[H])C([H])([H])C([H])([H])[H]':   'I'
}

D_CtermfreeAA_smarts = {
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])(C([H])([H])[H])O[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'N',
'[H]N([H])C(=O)C([H])([H])C([H])([H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'Q',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])(C([H])([H])[H])C(=O)O[H]':   'A',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C(=O)O[H])C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)N([H])C([H])(C([S-])([H])[H])C(=O)O[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])([H])C(=O)O[H]':   'G',
'[H]C([H])([H])C(=O)N1C([H])([H])C([H])([H])C([H])([H])C1([H])C(=O)O[H]':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C(=O)O[H])N([H])C(=O)C([H])([H])[H]':   'I',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(=N[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])[N+]([H])=C(N([H])[H])N([H])[H]':   'R',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])n([H])c(n1)[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])nc([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])[n+]([H])c([H])n1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])n([H])c([H])[n+]1[H]':   'H',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H]':   'K',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C(=O)O[H]':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])C([O-])=O':   'E',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C(=O)O[H]':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([O-])=O':   'D',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])O[H]':   'S',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])(C([H])([H])[H])O[H]':   'T',
'[H]N([H])C(=O)C([H])([H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'N',
'[H]N([H])C(=O)C([H])([H])C([H])([H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'Q',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])S[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])[H]':   'A',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'V',
'[H]C([H])([H])C([H])(C([H])([H])[H])C([H])([H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'L',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])C([H])([H])SC([H])([H])[H]':   'M',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])c([H])c([H])c([H])c1[H]':   'F',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])c([H])c([O-])c([H])c1[H]':   'Y',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([H])([H])c1c([H])n([H])c(c12)c([H])c([H])c([H])c2[H]':   'W',
'[H]C([H])([H])C(=O)N([H])C([H])(C([O-])=O)C([S-])([H])[H]':   'C',
'[H]C([H])([H])C(=O)N([H])C([H])([H])C([O-])=O':   'G',
'[H]C([H])([H])C(=O)N1C([H])([H])C([H])([H])C([H])([H])C1([H])C([O-])=O':   'P',
'[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])C([H])(C([O-])=O)N([H])C(=O)C([H])([H])[H]':   'I'
}
