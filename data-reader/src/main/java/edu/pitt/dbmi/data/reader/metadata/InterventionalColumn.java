/*
 * Copyright (C) 2018 University of Pittsburgh.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package edu.pitt.dbmi.data.reader.metadata;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * Dec 20, 2018 11:42:01 AM
 *
 * @author Kevin V. Bui (kvb2@pitt.edu)
 */
public class InterventionalColumn {

    @JsonProperty("value")
    private ColumnMetadata valueColumn;

    @JsonProperty("status")
    private ColumnMetadata statusColumn;

    public InterventionalColumn() {
    }

    public InterventionalColumn(ColumnMetadata valueColumn, ColumnMetadata statusColumn) {
        this.valueColumn = valueColumn;
        this.statusColumn = statusColumn;
    }

    @Override
    public String toString() {
        return "InterventionalColumn{" + "valueColumn=" + this.valueColumn + ", statusColumn=" + this.statusColumn + '}';
    }

    public ColumnMetadata getValueColumn() {
        return this.valueColumn;
    }

    public void setValueColumn(ColumnMetadata valueColumn) {
        this.valueColumn = valueColumn;
    }

    public ColumnMetadata getStatusColumn() {
        return this.statusColumn;
    }

    public void setStatusColumn(ColumnMetadata statusColumn) {
        this.statusColumn = statusColumn;
    }

}
