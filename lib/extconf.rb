require 'mkmf'

#* extconf.rb
#* Copyright (C) 2010  Embian Inc.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or 
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with This program.  If not, see <http://www.gnu.org/licenses/>.

EXT_COMPILE_C   = "#{COMPILE_C} -o $@"
EXT_COMPILE_CC  = "#{COMPILE_CXX} -o $@"
EXT_COMPILE_CXX = "#{COMPILE_CXX} -o $@"
EXT_COMPILE_CPP = "#{COMPILE_CXX} -o $@"

EXT_COMPILE_MAP = {'c'=>EXT_COMPILE_C, 'cc'=>EXT_COMPILE_CC, 'cxx'=>EXT_COMPILE_CXX, 'cpp'=>EXT_COMPILE_CPP}
def _configuration
  boostdir = "#{_dir_config('boost')}"
  etirmdir = "#{_dir_config('etirm')}"
  scppntdir = "#{_dir_config('scppnt')}"
  uncmindir = "#{_dir_config('uncmin')}"

  _error = "Required --with-boost-dir={boost_include}\n" unless File.exist? boostdir
  _error = "Required --with-etirm-dir={etirm_include}\n" unless File.exist? etirmdir
  _error = "Required --with-scppnt-dir={scppnt_include}\n" unless File.exist? scppntdir
  _error = "Required --with-uncmin-dir={uncmindir_include}\n" unless File.exist? uncmindir

  if _error
    message _error
    return false
  end
  
  boostdir = "#{boostdir}/"
  etirmdir = "#{etirmdir}/"
  scppntdir = "#{scppntdir}/"
  uncmindir = "#{uncmindir}/"
  
  puts "Creating SWIG wrapping modules"
  #   puts "swig -c++ -ruby -I. -I../include/etirm/src/ -I/usr/lib/ruby/1.8/x86_64-linux -o ./irt_ruby_warp.cpp ./irt_ruby.i"
  puts "swig -c++ -ruby -I. -I#{etirmdir} -I#{($extmk ? CONFIG["topdir"] : $topdir).quote} -o ./irt_ruby_warp.cpp ./irt_ruby.i"
  system "swig -c++ -ruby -I. -I#{etirmdir} -I#{($extmk ? CONFIG["topdir"] : $topdir).quote} -o ./irt_ruby_warp.cpp ./irt_ruby.i"
  
  CONFIG['CC'] = 'g++'
  CONFIG['LDSHARED'] = CONFIG['LDSHARED'].gsub(/^cc/, "$(CC)")
  
  CONFIG['boostdir'] = boostdir
  CONFIG['etirmdir'] = etirmdir
  CONFIG['scppntdir'] = scppntdir
  CONFIG['uncmindir'] = uncmindir

  dir_map = {}
  dir_map[$srcdir] = ""
  dir_map[boostdir] = "$(boostdir)"
  dir_map[etirmdir] = "$(etirmdir)"
  dir_map[scppntdir] = "$(scppntdir)"
  dir_map[uncmindir] = "$(uncmindir)"

  CONFIG['DEFS'] = "-DSCPPNT_BOUNDS_CHECK -DSCPPNT_NO_DIR_PREFIX -DETIRM_NO_DIR_PREFIX #{CONFIG['DEFS']}"

  $CFLAGS = "-Wno-deprecated #{CONFIG['CFLAGS']}"
  $CPPFLAGS = "-I$(boostdir) -I$(etirmdir) -I$(scppntdir) -I$(uncmindir) #{$CPPFLAGS}"
  
  $objs = []
  $srcs = []
  
  [$srcdir, etirmdir, scppntdir, uncmindir].each do |dir|
    srcs = Dir[File.join(dir, "*.{#{SRC_EXT.join(%q{,})}}")]
    prefix = dir_map[dir]
    srcs.collect!{|f| "#{prefix}#{File.basename(f)}"}
    for f in srcs
      obj = File.basename(f, ".*") << ".o"
      $objs.push(obj) unless $objs.index(obj)
    end
    $srcs.concat(srcs)
  end
  
  return true
end

# def _create_makefile(target, srcprefix = nil)
#   $target = target
#   libpath = $DEFLIBPATH|$LIBPATH
#   message "creating Makefile\n"
#   rm_f "conftest*"
#   if CONFIG["DLEXT"] == $OBJEXT
#     for lib in libs = $libs.split
#       lib.sub!(/-l(.*)/, %%"lib\\1.#{$LIBEXT}"%)
#     end
#     $defs.push(format("-DEXTLIB='%s'", libs.join(",")))
#   end
# 
#   if target.include?('/')
#     target_prefix, target = File.split(target)
#     target_prefix[0,0] = '/'
#   else
#     target_prefix = ""
#   end
# 
#   srcprefix ||= '$(srcdir)'
#   Config::expand(srcdir = srcprefix.dup)
# 
#   for i in $objs
#     i.sub!(/\.o\z/, ".#{$OBJEXT}")
#   end
# 
#   target = nil if $objs == ""
# 
#   if target and EXPORT_PREFIX
#     if File.exist?(File.join(srcdir, target + '.def'))
#       deffile = "$(srcdir)/$(TARGET).def"
#       unless EXPORT_PREFIX.empty?
#         makedef = %{-pe "sub!(/^(?=\\w)/,'#{EXPORT_PREFIX}') unless 1../^EXPORTS$/i"}
#       end
#     else
#       makedef = %{-e "puts 'EXPORTS', '#{EXPORT_PREFIX}Init_$(TARGET)'"}
#     end
#     if makedef
#       $distcleanfiles << '$(DEFFILE)'
#       origdef = deffile
#       deffile = "$(TARGET)-$(arch).def"
#     end
#   end
#   origdef ||= ''
# 
#   libpath = libpathflag(libpath)
# 
#   dllib = target ? "$(TARGET).#{CONFIG['DLEXT']}" : ""
#   staticlib = target ? "$(TARGET).#$LIBEXT" : ""
#   mfile = open("Makefile", "wb")
#   mfile.print configuration(srcprefix)
#   mfile.print "
# libpath = #{($DEFLIBPATH|$LIBPATH).join(" ")}
# LIBPATH = #{libpath}
# DEFFILE = #{deffile}
# 
# CLEANFILES = #{$cleanfiles.join(' ')}
# DISTCLEANFILES = #{$distcleanfiles.join(' ')}
# 
# extout = #{$extout}
# extout_prefix = #{$extout_prefix}
# target_prefix = #{target_prefix}
# LOCAL_LIBS = #{$LOCAL_LIBS}
# LIBS = #{$LIBRUBYARG} #{$libs} #{$LIBS}
# SRCS = #{$srcs.join(" ")}
# OBJS = #{$objs.join(" ")}
# TARGET = #{target}
# DLLIB = #{dllib}
# EXTSTATIC = #{$static || ""}
# STATIC_LIB = #{staticlib unless $static.nil?}
# #{!$extout && defined?($installed_list) ? "INSTALLED_LIST = #{$installed_list}\n" : ""}
# "
#   install_dirs.each {|d| mfile.print("%-14s= %s\n" % d) if /^[[:upper:]]/ =~ d[0]}
#   n = ($extout ? '$(RUBYARCHDIR)/' : '') + '$(TARGET).'
#   mfile.print "
# TARGET_SO     = #{($extout ? '$(RUBYARCHDIR)/' : '')}$(DLLIB)
# CLEANLIBS     = #{n}#{CONFIG['DLEXT']} #{n}il? #{n}tds #{n}map
# CLEANOBJS     = *.#{$OBJEXT} *.#{$LIBEXT} *.s[ol] *.pdb *.exp *.bak #{$objs.join(" ")}
# 
# all:            #{$extout ? "install" : target ? "$(DLLIB)" : "Makefile"}
# static:         $(STATIC_LIB)#{$extout ? " install-rb" : ""}
# "
#   mfile.print CLEANINGS
#   dirs = []
#   mfile.print "install: install-so install-rb\n\n"
#   sodir = (dir = "$(RUBYARCHDIR)").dup
#   mfile.print("install-so: #{dir}\n")
#   if target
#     f = "$(DLLIB)"
#     dest = "#{dir}/#{f}"
#     mfile.print "install-so: #{dest}\n"
#     unless $extout
#       mfile.print "#{dest}: #{f}\n"
#       if (sep = config_string('BUILD_FILE_SEPARATOR'))
#         f.gsub!("/", sep)
#         dir.gsub!("/", sep)
#         sep = ":/="+sep
#         f.gsub!(/(\$\(\w+)(\))/) {$1+sep+$2}
#         f.gsub!(/(\$\{\w+)(\})/) {$1+sep+$2}
#         dir.gsub!(/(\$\(\w+)(\))/) {$1+sep+$2}
#         dir.gsub!(/(\$\{\w+)(\})/) {$1+sep+$2}
#       end
#       mfile.print "\t$(INSTALL_PROG) #{f} #{dir}\n"
#       if defined?($installed_list)
#         mfile.print "\t@echo #{dir}/#{File.basename(f)}>>$(INSTALLED_LIST)\n"
#       end
#     end
#   end
#   mfile.print("install-rb: pre-install-rb install-rb-default\n")
#   mfile.print("install-rb-default: pre-install-rb-default\n")
#   mfile.print("pre-install-rb: Makefile\n")
#   mfile.print("pre-install-rb-default: Makefile\n")
#   for sfx, i in [["-default", [["lib/**/*.rb", "$(RUBYLIBDIR)", "lib"]]], ["", $INSTALLFILES]]
#     files = install_files(mfile, i, nil, srcprefix) or next
#     for dir, *files in files
#       unless dirs.include?(dir)
#         dirs << dir
#         mfile.print "pre-install-rb#{sfx}: #{dir}\n"
#       end
#       files.each do |f|
#         dest = "#{dir}/#{File.basename(f)}"
#         mfile.print("install-rb#{sfx}: #{dest}\n")
#         mfile.print("#{dest}: #{f}\n\t$(#{$extout ? 'COPY' : 'INSTALL_DATA'}) ")
#         sep = config_string('BUILD_FILE_SEPARATOR')
#         if sep
#           f = f.gsub("/", sep)
#           sep = ":/="+sep
#           f = f.gsub(/(\$\(\w+)(\))/) {$1+sep+$2}
#           f = f.gsub(/(\$\{\w+)(\})/) {$1+sep+$2}
#         else
#           sep = ""
#         end
#         mfile.print("#{f} $(@D#{sep})\n")
#         if defined?($installed_list) and !$extout
#           mfile.print("\t@echo #{dest}>>$(INSTALLED_LIST)\n")
#         end
#       end
#     end
#   end
#   dirs.unshift(sodir) if target and !dirs.include?(sodir)
#   dirs.each {|dir| mfile.print "#{dir}:\n\t$(MAKEDIRS) $@\n"}
# 
#   mfile.print "\nsite-install: site-install-so site-install-rb\nsite-install-so: install-so\nsite-install-rb: install-rb\n\n"
# 
#   return unless target
# 
#   mfile.puts SRC_EXT.collect {|ext| ".path.#{ext} = $(VPATH)"} if $nmake == ?b
#   mfile.print ".SUFFIXES: .#{SRC_EXT.join(' .')} .#{$OBJEXT}\n"
#   mfile.print "\n"
#   exts = [].concat(CXX_EXT).concat(%w[c])
#   exts.each do |ext|
#     ext_compile_rule = EXT_COMPILE_MAP[ext.downcase]
#     raise "Can't find compile rule of '#{ext}' source file" unless ext_compile_rule
#     COMPILE_RULES.each do |rule|
#       mfile.printf(rule, ext, $OBJEXT)
#       mfile.printf("\n\t%s\n\n", ext_compile_rule)
#     end
#   end
# 
#   mfile.print "$(RUBYARCHDIR)/" if $extout
#   mfile.print "$(DLLIB): "
#   mfile.print "$(DEFFILE) " if makedef
#   mfile.print "$(OBJS) Makefile\n"
#   mfile.print "\t@-$(RM) $@\n"
#   mfile.print "\t@-$(MAKEDIRS) $(@D)\n" if $extout
#   link_so = LINK_SO.gsub(/^/, "\t")
#   mfile.print link_so, "\n\n"
#   unless $static.nil?
#     mfile.print "$(STATIC_LIB): $(OBJS)\n\t"
#     mfile.print "$(AR) #{config_string('ARFLAGS') || 'cru '}$@ $(OBJS)"
#     config_string('RANLIB') do |ranlib|
#       mfile.print "\n\t@-#{ranlib} $(DLLIB) 2> /dev/null || true"
#     end
#   end
#   mfile.print "\n\n"
#   if makedef
#     mfile.print "$(DEFFILE): #{origdef}\n"
#     mfile.print "\t$(RUBY) #{makedef} #{origdef} > $@\n\n"
#   end
# 
#   headers = %w[ruby.h defines.h]
#   if RULE_SUBST
#     headers.each {|h| h.sub!(/.*/) {|*m| RULE_SUBST % m}}
#   end
#   headers << $config_h if $config_h
#   headers << "$(RUBY_EXTCONF_H)" if $extconf_h
#   mfile.print "$(OBJS): ", headers.join(' '), "\n"
# 
#   $makefile_created = true
# ensure
#   mfile.close if mfile
# end

def _dir_config(target, idefault=nil, ldefault=nil)
  if dir = with_config(target + "-dir", (idefault unless ldefault))
    defaults = Array === dir ? dir : dir.split(File::PATH_SEPARATOR)
    idefault = ldefault = nil
  end
  return dir
end

# _create_makefile("irt_ruby") if _configuration
create_makefile("irt_ruby") if _configuration
