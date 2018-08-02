#ifndef RESULT_VIEWER_H
#define RESULT_VIEWER_H

#include <QtWidgets>
#include "integrationresult.h"

namespace Integration {
    template <typename T>
    void view_result(const Result<T>& result) {
        int argc = 1;
        char *argv[] = { const_cast<char *>("Result")};
        QApplication app(argc, argv);

        const unsigned int MAXSIZE = 700;
        const auto& [xlow, xhigh] = result.origin->intervals[0];
        const auto& [ylow, yhigh] = result.origin->intervals[1];
        auto xmid = (xhigh + xlow) / 2, ymid = (yhigh  + ylow) / 2;
        auto scale = MAXSIZE / std::max(xhigh - xlow, yhigh - ylow);

        auto pix = new QPixmap((xhigh - xlow) * scale, (yhigh - ylow) * scale);
        pix->fill(Qt::white);

        auto painter = new QPainter(pix);
        painter->setPen(QColor(0, 0, 0, 200));
        painter->translate(pix->width() / 2, pix->height() / 2);
        painter->scale(1, -1);

        auto label = new QLabel;
        label->setFixedSize(pix->size());

        QBrush border_brush(QColor("#2192d3"));
        QBrush interior_brush(QColor("#7bc1e9"));
        QBrush no_brush = Qt::NoBrush;

        for (const auto& [hc, state] : result.cubes) {
            QPainterPath path;
            QBrush brush = no_brush;

            if (state == CONTAINED)
                brush = interior_brush;
            else if (state == INDEFINITE)
                brush = border_brush;

            const auto& [xlow, xhigh] = hc->intervals[0];
            const auto& [ylow, yhigh] = hc->intervals[1];

            auto cube = QRectF((xlow - xmid) * scale, (ylow - ymid) * scale,
                               (xhigh - xlow) * scale, (yhigh - ylow) * scale);

            path.addRect(cube);
            painter->fillPath(path, brush);
            painter->drawPath(path);
        }

        label->setPixmap(*pix);
        label->show();
        app.exec();
    }
}

#endif // RESULT_VIEWER_H
