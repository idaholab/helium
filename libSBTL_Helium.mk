#
# build libSBTL_Helium using libMesh's build system
#
LIBSBTL_HELIUM_DIR       := $(HELIUM_DIR)/contrib/libSBTL_Helium

LIBSBTL_HELIUM_srcfiles  :=
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/CP_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/CV_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/ETA_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/G_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/HS_FLASH.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/HV_FLASH.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/LAMBDA_VU_HE.cpp
#LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/LibSBTL_vu_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/P_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/PH_FLASH.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/PS_FLASH.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/PT_FLASH.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/S_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/T_VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/U_VH_HE_INI.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/U_VP_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/VU_HP_HE_INI.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/VU_HE.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/VU_SH_HE_INI.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/VU_SP_HE_INI.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/VU_TP_HE_INI.cpp
LIBSBTL_HELIUM_srcfiles  += $(LIBSBTL_HELIUM_DIR)/W_VU_HE.cpp

LIBSBTL_HELIUM_objects   := $(patsubst %.cpp, %.$(obj-suffix), $(LIBSBTL_HELIUM_srcfiles))
LIBSBTL_HELIUM_deps      := $(patsubst %.$(obj-suffix), %.$(obj-suffix).d, $(LIBSBTL_HELIUM_objects))
LIBSBTL_HELIUM_LIB       := $(LIBSBTL_HELIUM_DIR)/libSBTL_Helium-$(METHOD).la

app_INCLUDES += -I$(LIBSBTL_HELIUM_DIR)
app_LIBS += $(LIBSBTL_HELIUM_LIB)

$(LIBSBTL_HELIUM_LIB): $(LIBSBTL_HELIUM_objects)
	@echo "Linking Library "$@"..."
	@$(libmesh_LIBTOOL) --tag=CC $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CC) $(libmesh_CFLAGS) -o $@ $(LIBSBTL_HELIUM_objects) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(LIBSBTL_HELIUM_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(LIBSBTL_HELIUM_LIB) $(LIBSBTL_HELIUM_DIR)

$(app_EXEC): $(LIBSBTL_HELIUM_LIB)

-include $(LIBSBTL_HELIUM_deps)

cleanlibsbtl_nitrogen:
	@echo "Cleaning libSBTL_Helium"
	@rm -f $(LIBSBTL_HELIUM_objects)
	@rm -f $(LIBSBTL_HELIUM_deps)
	@rm -f $(LIBSBTL_HELIUM_LIB)
	@rm -f $(LIBSBTL_HELIUM_DIR)/libSBTL_Helium-$(METHOD)*.dylib
	@rm -f $(LIBSBTL_HELIUM_DIR)/libSBTL_Helium-$(METHOD)*.so*
	@rm -f $(LIBSBTL_HELIUM_DIR)/libSBTL_Helium-$(METHOD)*.a
