// Module:  Log4CPLUS
// File:    factory.cxx
// Created: 2/2002
// Author:  Tad E. Smith
//
//
// Copyright 2002-2009 Tad E. Smith
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "dcmtk/oflog/spi/factory.h"
#include "dcmtk/oflog/spi/logfact.h"
#include "dcmtk/oflog/consap.h"
#include "dcmtk/oflog/fileap.h"
#include "dcmtk/oflog/nullap.h"
#include "dcmtk/oflog/socketap.h"
#include "dcmtk/oflog/syslogap.h"
#include "dcmtk/oflog/helpers/loglog.h"
#include "dcmtk/oflog/helpers/threads.h"

#if defined (_WIN32)
#  if defined (LOG4CPLUS_HAVE_NT_EVENT_LOG)
#    include "dcmtk/oflog/ntelogap.h"
#  endif
#  if defined (LOG4CPLUS_HAVE_WIN32_CONSOLE)
#    include "dcmtk/oflog/winconap.h"
#  endif
#  include "dcmtk/oflog/windebap.h"
#endif


using namespace log4cplus;
using namespace log4cplus::helpers;
using namespace log4cplus::spi;


///////////////////////////////////////////////////////////////////////////////
// LOCAL file class definitions
///////////////////////////////////////////////////////////////////////////////

namespace log4cplus {

namespace spi {

    BaseFactory::~BaseFactory()
    { }


    AppenderFactory::AppenderFactory()
    { }

    AppenderFactory::~AppenderFactory()
    { }


    LayoutFactory::LayoutFactory()
    { }

    LayoutFactory::~LayoutFactory()
    { }


    FilterFactory::FilterFactory()
    { }

    FilterFactory::~FilterFactory()
    { }


    LoggerFactory::~LoggerFactory()
    { }


} // namespace spi


namespace
{


template <typename ProductFactoryBase>
class LocalFactoryBase
    : public ProductFactoryBase
{
public:
    LocalFactoryBase (tchar const * n)
        : name (n)
    { }

    virtual log4cplus::tstring getTypeName()
    {
        return name;
    }

private:
    log4cplus::tstring name;
};


template <typename LocalProduct, typename ProductFactoryBase>
class FactoryTempl
    : public LocalFactoryBase<ProductFactoryBase>
{
public:
    typedef typename ProductFactoryBase::ProductPtr ProductPtr;

    FactoryTempl (tchar const * n)
        : LocalFactoryBase<ProductFactoryBase> (n)
    { }

    virtual ProductPtr createObject (Properties const & props, log4cplus::tstring& error)
    {
        error.clear();
        return ProductPtr (new LocalProduct (props, error));
    }
};


} // namespace


#define REG_PRODUCT(reg, productprefix, productname, productns, productfact) \
reg.put (                                                               \
    OFauto_ptr<productfact> (                                        \
        new FactoryTempl<productns productname, productfact> (          \
            LOG4CPLUS_TEXT(productprefix)                               \
            LOG4CPLUS_TEXT(#productname))))


#define REG_APPENDER(reg, appendername)                             \
REG_PRODUCT (reg, "log4cplus::", appendername, log4cplus::, AppenderFactory)

#define REG_LAYOUT(reg, layoutname)                                 \
REG_PRODUCT (reg, "log4cplus::", layoutname, log4cplus::, LayoutFactory)

#define REG_FILTER(reg, filtername)                                 \
REG_PRODUCT (reg, "log4cplus::spi::", filtername, spi::, FilterFactory)


void initializeFactoryRegistry()
{
    AppenderFactoryRegistry& reg = getAppenderFactoryRegistry();
    REG_APPENDER (reg, ConsoleAppender);
    REG_APPENDER (reg, NullAppender);
    REG_APPENDER (reg, FileAppender);
    REG_APPENDER (reg, RollingFileAppender);
    REG_APPENDER (reg, DailyRollingFileAppender);
    REG_APPENDER (reg, SocketAppender);
#if defined(_WIN32) && !defined(__MINGW32__)
#if defined(LOG4CPLUS_HAVE_NT_EVENT_LOG)
    REG_APPENDER (reg, NTEventLogAppender);
#  endif
#  if defined(LOG4CPLUS_HAVE_WIN32_CONSOLE)
    REG_APPENDER (reg, Win32ConsoleAppender);
#  endif
    REG_APPENDER (reg, Win32DebugAppender);
#elif defined(LOG4CPLUS_HAVE_SYSLOG_H)
    REG_APPENDER (reg, SysLogAppender);
#endif // defined(_WIN32) && !defined(__MINGW32__)

    LayoutFactoryRegistry& reg2 = getLayoutFactoryRegistry();
    REG_LAYOUT (reg2, SimpleLayout);
    REG_LAYOUT (reg2, TTCCLayout);
    REG_LAYOUT (reg2, PatternLayout);

    FilterFactoryRegistry& reg3 = getFilterFactoryRegistry();
    REG_FILTER (reg3, DenyAllFilter);
    REG_FILTER (reg3, LogLevelMatchFilter);
    REG_FILTER (reg3, LogLevelRangeFilter);
    REG_FILTER (reg3, StringMatchFilter);
}




///////////////////////////////////////////////////////////////////////////////
// public methods
///////////////////////////////////////////////////////////////////////////////

namespace spi
{


AppenderFactoryRegistry&
getAppenderFactoryRegistry()
{
    static AppenderFactoryRegistry singleton;
    return singleton;
}


LayoutFactoryRegistry&
getLayoutFactoryRegistry()
{
    static LayoutFactoryRegistry singleton;
    return singleton;
}


FilterFactoryRegistry&
getFilterFactoryRegistry()
{
    static FilterFactoryRegistry singleton;
    return singleton;
}


} // namespace spi


} // namespace log4cplus
