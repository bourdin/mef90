diff -r 55e134b824e3 config/packages/exodusii.py
--- a/config/packages/exodusii.py	Thu Nov 29 12:01:56 2012 -0600
+++ b/config/packages/exodusii.py	Sat Nov 16 19:25:43 2013 -0600
@@ -6,7 +6,7 @@
 class Configure(config.package.Package):
   def __init__(self, framework):
     config.package.Package.__init__(self, framework)
-    self.download   = ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/exodusii-5.14-petsc.tgz']
+    self.download   = ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/exodusii-5.22b.tar.gz']
     self.liblist    = [['libexoIIv2for.a', 'libexodus.a'], ['libexoIIv2for.a', 'libexoIIv2c.a']]
     self.functions  = ['ex_close']
     self.includes   = ['exodusII.h']
@@ -21,30 +21,48 @@
     return
 
   def Install(self):
-    self.framework.log.write('exodusIIDir = '+self.packageDir+' installDir '+self.installDir+'\n')
+    self.logPrintBox('Compiling ExodusII; this may take several minutes')
+    import os
+    import shutil
+    configOpts     = []
+    configOpts.append('RANLIB="'+self.setCompilers.RANLIB+'"')
+    configOpts.append('AR="'+self.setCompilers.AR+' '+self.setCompilers.AR_FLAGS+'"')
+    configOpts.append('NETCDF="'+self.installDir+'"')
+
+    self.setCompilers.pushLanguage('C')
+    configOpts.append('CC="'+self.setCompilers.getCompiler()+'"')
+    configOpts.append('CCOPTIONS="'+self.setCompilers.getCompilerFlags()+' -DADDC_ "')
+    self.setCompilers.popLanguage()
+
+    if hasattr(self.setCompilers, 'FC'):
+      self.setCompilers.pushLanguage('FC')
+      configOpts.append('FC="'+self.setCompilers.getCompiler()+'"')
+      configOpts.append('F77OPTIONS="'+self.setCompilers.getCompilerFlags()+'"')
+    self.setCompilers.popLanguage()
 
     mkfile = 'make.inc'
     g = open(os.path.join(self.packageDir, mkfile), 'w')
     self.framework.log.write(repr(dir(self.setCompilers)))
-    self.setCompilers.pushLanguage('C')
-    g.write('CC = '+self.setCompilers.getCompiler()+'\n')
-    g.write('CC_FLAGS = '+self.setCompilers.getCompilerFlags()+'\n')
-    self.setCompilers.popLanguage()
 
-    self.setCompilers.pushLanguage('FC')
-    g.write('FC = '+self.setCompilers.getCompiler()+'\n')
-    g.write('FC_FLAGS = '+self.setCompilers.getCompilerFlags()+'\n')
-    self.setCompilers.popLanguage()
-    g.write('RANLIB      = '+self.setCompilers.RANLIB+'\n')
-    g.write('AR      = '+self.setCompilers.AR+'\n')
-    g.write('AR_FLAGS      = '+self.setCompilers.AR_FLAGS+'\n')
+    args = ' '.join(configOpts)
+    fd = file(os.path.join(self.packageDir,'exodusii'), 'w')
+    fd.write(args)
+    fd.close()
 
-    g.close()
-
-    if self.installNeeded(mkfile):
+    if self.installNeeded('exodusii'):
+      cincludes  = ['exodusII.h','exodusII_cfg.h','exodusII_int.h','exodusII_par.h']
+      fincludes  = ['exodusII.inc','exodusII_int.inc']
       try:
         self.logPrintBox('Compiling ExodusII; this may take several minutes')
-        output,err,ret = config.base.Configure.executeShellCommand('cd '+self.packageDir+' && make -f Makefile.petsc clean && make -f Makefile.petsc && make -f Makefile.petsc install', timeout=2500, log = self.framework.log)
+        output,err,ret = config.base.Configure.executeShellCommand('cd '+self.packageDir+' && make -f Makefile.standalone libexodus.a '+args, timeout=2500, log = self.framework.log)
+        shutil.copy(os.path.join(self.packageDir,'libexodus.a'),os.path.join(self.installDir,'lib'))
+        for i in cincludes:
+          shutil.copy(os.path.join(self.packageDir,'cbind','include',i),os.path.join(self.installDir,'include'))
+        if hasattr(self.setCompilers, 'FC'):
+          output,err,ret = config.base.Configure.executeShellCommand('cd '+self.packageDir+' && make -f Makefile.standalone libexoIIv2for.a '+args, timeout=2500, log = self.framework.log)
+          shutil.copy(os.path.join(self.packageDir,'libexoIIv2for.a'),os.path.join(self.installDir,'lib'))
+          for i in fincludes:
+            shutil.copy(os.path.join(self.packageDir,'forbind','include',i),os.path.join(self.installDir,'include'))
       except RuntimeError, e:
         raise RuntimeError('Error running make on exodusII: '+str(e))
       self.postInstall(output+err, mkfile)
