diff -r a3b0a0b09619 config/packages/netcdf.py
--- a/config/packages/netcdf.py	Wed Jun 13 16:26:42 2012 -0500
+++ b/config/packages/netcdf.py	Tue Jul 31 21:41:37 2012 -0500
@@ -8,18 +8,20 @@ class Configure(config.package.GNUPackag
     config.package.GNUPackage.__init__(self, framework)
     self.downloadpath    = 'http://www.unidata.ucar.edu/downloads/netcdf/ftp/'
     self.downloadext     = 'tar.gz'
-    self.downloadversion = '4.1.1'
+    self.downloadversion = '4.2.1'
     self.functions       = ['nccreate']
     self.includes        = ['netcdf.h']
-    self.liblist         = [['libnetcdf_c++.a','libnetcdf.a']]
+    #self.liblist         = [['libnetcdf_c++.a','libnetcdf.a']]
+    self.liblist         = [['libnetcdf.a']]
     self.cxx             = 1
     return
 
   def setupDependencies(self, framework):
     config.package.GNUPackage.setupDependencies(self, framework)
     self.mpi   = framework.require('config.packages.MPI', self)
-    self.hdf5  = framework.require('config.packages.hdf5', self)
-    self.odeps = [self.mpi, self.hdf5]
+    #self.hdf5  = framework.require('config.packages.hdf5', self)
+    #self.odeps = [self.mpi, self.hdf5]
+    self.odeps = [self.mpi]
     return
 
   def Install(self):
@@ -41,6 +43,7 @@ class Configure(config.package.GNUPackag
     configOpts.append('--prefix='+self.installDir)
     configOpts.append('--libdir='+os.path.join(self.installDir,self.libdir))
     configOpts.append('--disable-dap')
+    configOpts.append('--disable-netcdf-4')
 
     self.setCompilers.pushLanguage('C')
     cflags = self.setCompilers.getCompilerFlags().replace('-Wall','').replace('-Wshadow','')
