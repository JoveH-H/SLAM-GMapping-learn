#ifndef ACCESSTATE_H
#define ACCESSTATE_H

namespace GMapping
{
    /* 可访问性状态 枚举类型的定义 */
    enum AccessibilityState
    {
        Outside = 0x0,  /* 外部 */
        Inside = 0x1,   /* 内部 */
        Allocated = 0x2 /* 分配 */
    };
};

#endif
