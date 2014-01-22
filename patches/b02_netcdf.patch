diff -r 55e134b824e3 config/packages/netcdf.py
--- a/config/packages/netcdf.py	Thu Nov 29 12:01:56 2012 -0600
+++ b/config/packages/netcdf.py	Sat Nov 16 19:26:07 2013 -0600
@@ -8,10 +8,10 @@
     config.package.GNUPackage.__init__(self, framework)
     self.downloadpath    = 'http://www.unidata.ucar.edu/downloads/netcdf/ftp/'
     self.downloadext     = 'tar.gz'
-    self.downloadversion = '4.1.1'
+    self.downloadversion = '4.2.1.1'
     self.functions       = ['nccreate']
     self.includes        = ['netcdf.h']
-    self.liblist         = [['libnetcdf_c++.a','libnetcdf.a']]
+    self.liblist         = [['libnetcdf.a']]
     self.cxx             = 1
     return
 
@@ -20,45 +20,36 @@
     self.mpi   = framework.require('config.packages.MPI', self)
     self.hdf5  = framework.require('config.packages.hdf5', self)
     self.odeps = [self.mpi, self.hdf5]
+    self.deps  = [self.mpi]
     return
 
   def Install(self):
     import os, sys
 
-    makeinc        = os.path.join(self.packageDir, 'make.inc')
-    installmakeinc = os.path.join(self.confDir, 'NetCDF')
-    configEnv      = []
     configOpts     = []
     # Unused flags: F90, CPPFLAGS, LIBS, FLIBS
-    g = open(makeinc, 'w')
-    g.write('AR             = '+self.setCompilers.AR+'\n')
-    g.write('ARFLAGS        = '+self.setCompilers.AR_FLAGS+'\n')
-    configEnv.append('AR="'+self.setCompilers.AR+'"')
-    configEnv.append('ARFLAGS="'+self.setCompilers.AR_FLAGS+'"')
+    configOpts.append('AR="'+self.setCompilers.AR+'"')
+    configOpts.append('ARFLAGS="'+self.setCompilers.AR_FLAGS+'"')
 
-    g.write('NETCDF_ROOT    = '+self.packageDir+'\n')
-    g.write('PREFIX         = '+self.installDir+'\n')
     configOpts.append('--prefix='+self.installDir)
     configOpts.append('--libdir='+os.path.join(self.installDir,self.libdir))
     configOpts.append('--disable-dap')
+    configOpts.append('--disable-hdf4')
+    configOpts.append('--disable-netcdf-4')
 
     self.setCompilers.pushLanguage('C')
     cflags = self.setCompilers.getCompilerFlags().replace('-Wall','').replace('-Wshadow','')
     cflags += ' ' + self.headers.toString(self.mpi.include)+' '+self.headers.toString('.')
-    g.write('CC             = '+self.setCompilers.getCompiler()+'\n')
-    g.write('CFLAGS         = '+cflags+'\n')
-    configEnv.append('CC="'+self.setCompilers.getCompiler()+'"')
-    configEnv.append('CFLAGS="'+cflags+'"')
+    configOpts.append('CC="'+self.setCompilers.getCompiler()+'"')
+    configOpts.append('CFLAGS="'+cflags+'"')
     self.setCompilers.popLanguage()
 
     if hasattr(self.setCompilers, 'CXX'):
       self.setCompilers.pushLanguage('Cxx')
       cxxflags = self.setCompilers.getCompilerFlags().replace('-Wall','').replace('-Wshadow','')
       cxxflags += ' ' + self.headers.toString(self.mpi.include)+' '+self.headers.toString('.')
-      g.write('CXX            = '+self.setCompilers.getCompiler()+'\n')
-      g.write('CXXFLAGS       = '+cflags+'\n')
-      configEnv.append('CXX="'+self.setCompilers.getCompiler()+'"')
-      configEnv.append('CXXFLAGS="'+cxxflags+'"')
+      configOpts.append('CXX="'+self.setCompilers.getCompiler()+'"')
+      configOpts.append('CXXFLAGS="'+cxxflags+'"')
       self.setCompilers.popLanguage()
     else:
       configOpts.append('--disable-cxx')
@@ -67,12 +58,10 @@
       self.setCompilers.pushLanguage('FC')
       fcflags = self.setCompilers.getCompilerFlags().replace('-Wall','').replace('-Wshadow','')
       fcflags += ' ' + self.headers.toString(self.mpi.include)+' '+self.headers.toString('.')
-      g.write('FC             = '+self.setCompilers.getCompiler()+'\n')
-      g.write('FCFLAGS        = '+fcflags+'\n')
-      configEnv.append('FC="'+self.setCompilers.getCompiler()+'"')
-      configEnv.append('FCFLAGS="'+fcflags+'"')
+      configOpts.append('FC="'+self.setCompilers.getCompiler()+'"')
+      configOpts.append('FCFLAGS="'+fcflags+'"')
       if self.compilers.fortranIsF90:
-        configEnv.append('F90="'+self.setCompilers.getCompiler()+'"')
+        configOpts.append('F90="'+self.setCompilers.getCompiler()+'"')
       else:
         configOpts.append('--disable-f90')
       self.setCompilers.popLanguage()
@@ -81,15 +70,19 @@
 
     if self.setCompilers.sharedLibraries:
       configOpts.append('--enable-shared')
-    g.close()
 
-    if self.installNeeded('make.inc'):    # Now compile & install
+    args = ' '.join(configOpts)
+    fd = file(os.path.join(self.packageDir,'netcdf'), 'w')
+    fd.write(args)
+    fd.close()
+
+    if self.installNeeded('netcdf'):
       try:
         self.logPrintBox('Configuring NetCDF; this may take several minutes')
-        output,err,ret  = self.executeShellCommand('cd '+self.packageDir+' && '+' '.join(configEnv)+' ./configure '+' '.join(configOpts), timeout=2500, log = self.framework.log)
+        output,err,ret  = self.executeShellCommand('cd '+self.packageDir+' && ./configure '+args, timeout=2500, log = self.framework.log)
         self.logPrintBox('Compiling & installing NetCDF; this may take several minutes')
         output,err,ret  = self.executeShellCommand('cd '+self.packageDir+' && make clean && make && make install && make clean', timeout=2500, log = self.framework.log)
       except RuntimeError, e:
         raise RuntimeError('Error running make on NetCDF: '+str(e))
-      self.postInstall(output+err,'make.inc')
+      self.postInstall(output+err,'netcdf')
     return self.installDir
