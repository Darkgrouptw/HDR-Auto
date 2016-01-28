/****************************************************************************
** Meta object code from reading C++ file 'hdrauto.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../hdrauto.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'hdrauto.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_HDRAuto_t {
    QByteArrayData data[8];
    char stringdata0[89];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_HDRAuto_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_HDRAuto_t qt_meta_stringdata_HDRAuto = {
    {
QT_MOC_LITERAL(0, 0, 7), // "HDRAuto"
QT_MOC_LITERAL(1, 8, 13), // "OpenFileEvent"
QT_MOC_LITERAL(2, 22, 0), // ""
QT_MOC_LITERAL(3, 23, 10), // "StartEvent"
QT_MOC_LITERAL(4, 34, 12), // "SetMaskEvent"
QT_MOC_LITERAL(5, 47, 12), // "Text1Changed"
QT_MOC_LITERAL(6, 60, 12), // "Text2Changed"
QT_MOC_LITERAL(7, 73, 15) // "BarChangedEvent"

    },
    "HDRAuto\0OpenFileEvent\0\0StartEvent\0"
    "SetMaskEvent\0Text1Changed\0Text2Changed\0"
    "BarChangedEvent"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_HDRAuto[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    0,   46,    2, 0x08 /* Private */,
       5,    1,   47,    2, 0x08 /* Private */,
       6,    1,   50,    2, 0x08 /* Private */,
       7,    0,   53,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,

       0        // eod
};

void HDRAuto::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        HDRAuto *_t = static_cast<HDRAuto *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->OpenFileEvent(); break;
        case 1: _t->StartEvent(); break;
        case 2: _t->SetMaskEvent(); break;
        case 3: _t->Text1Changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->Text2Changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->BarChangedEvent(); break;
        default: ;
        }
    }
}

const QMetaObject HDRAuto::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_HDRAuto.data,
      qt_meta_data_HDRAuto,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *HDRAuto::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *HDRAuto::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_HDRAuto.stringdata0))
        return static_cast<void*>(const_cast< HDRAuto*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int HDRAuto::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
